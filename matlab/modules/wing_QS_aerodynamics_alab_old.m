function [L_R,L_L, D_R, D_L, M_R, M_L, ...
          F_rot_R, F_rot_L, M_rot_R, M_rot_L, ...
          alpha_R, alpha_L, U_alpha_R_dot, U_alpha_L_dot, U_R, U_L]...
          = wing_QS_aerodynamics_old(INSECT, W_R, W_L, W_R_dot, W_L_dot, varargin)
% compute forces and moment from QS model
% [L_R L_L ...] = wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot)
% computes the aerodynamic forces with small velocity assumption
%
% [L_R L_L ...] = wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L)
% computes the aerodynamic forces while considering the effects of x_dot
% and W on the flow around the wing, by integrating infinitesimal forces
% along the wing span

global e1 e2 e3
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

if nargin < 6
    %% small x_dot, Omega assumpton
    U_R = INSECT.tilde_r_2*INSECT.l*cross(W_R,e2);
    U_R_dot = INSECT.tilde_r_2*INSECT.l*cross(W_R_dot,e2);
    [L_R D_R M_R]=compute_LD(INSECT, U_R);
    [F_rot_R M_rot_R alpha_R U_alpha_R_dot]=compute_rotational_force(INSECT, U_R, U_R_dot);
    
    U_L = -INSECT.tilde_r_2*INSECT.l*cross(W_L,e2);
    U_L_dot = -INSECT.tilde_r_2*INSECT.l*cross(W_L_dot,e2);
    [L_L D_L M_L]=compute_LD(INSECT, U_L);
    M_L = -M_L;
    [F_rot_L M_rot_L alpha_L U_alpha_L_dot]=compute_rotational_force(INSECT, U_L, U_L_dot);
    M_rot_L = -M_rot_L;
    %display('here?');
    
else
    %display('here?2');
    %% integrate infinitesimal forces over the wing span
    x_dot=varargin{1};
    R=varargin{2};
    W=varargin{3};
    Q_R=varargin{4};
    Q_L=varargin{5};
    Euler_R=varargin{6};
    Euler_R_dot=varargin{7};
    Euler_L=varargin{8};
    Euler_L_dot=varargin{9};
    
    N_r=50;
    rs=linspace(0,INSECT.l,N_r);
    dr=rs(2)-rs(1);
    L_R = zeros(3,1);
    D_R = zeros(3,1);
    M_R = zeros(3,1);
    L_L = zeros(3,1);
    D_L = zeros(3,1);
    M_L = zeros(3,1);

    %% Translational forces and moments
    % 
    for i=1:N_r-1
        r=rs(i);
        c=polyval(INSECT.cr_poly,r*1e2)*1e-2;
        xLE = polyval(INSECT.cr_LE_poly,r*1e2)*1e-2;
        xCP = (xLE - 0.25 * c) * e1 ;
        xPivot = 0 * e1; % pitch rotation at the root
        
        U_R = (eye(3)-e2*e2')*Q_R'*(R'*x_dot+hat(W)*INSECT.mu_R) + r*hat(Q_R*W+W_R)*e2;
        U_R_norm = U_R / norm(U_R);
        [L D M]=compute_LD(INSECT, U_R);
        
        L_R = L_R + L/INSECT.S*c*dr;
        D_R = D_R + D/INSECT.S*c*dr;
        M_R = M_R + cross(r*e2,L)/INSECT.S*c*dr ;
        M_R = M_R + cross(r*e2,D)/INSECT.S*c*dr ;
        M_R = M_R - cross(L, (xCP - xPivot))/INSECT.S*c*dr ;
        
        U_L = (eye(3)-e2*e2')*Q_L'*(R'*x_dot+hat(W)*INSECT.mu_L) - r*hat(Q_L*W+W_L)*e2;
        [L D M]=compute_LD(INSECT, U_L);
        
        L_L = L_L + L/INSECT.S*c*dr;
        D_L = D_L + D/INSECT.S*c*dr;
        M_L = M_L - cross(r*e2,(L))/INSECT.S*c*dr;
        M_L = M_L - cross(r*e2,(D))/INSECT.S*c*dr;
        M_L = M_L - cross(L, (xCP - xPivot))/INSECT.S*c*dr;
    end
    
    %% rotational forces and moments
    F_rot_R=zeros(3,1);
    M_rot_R=zeros(3,1);
    alpha_R=0;
    U_alpha_R_dot=0;
    
    F_rot_L=zeros(3,1);
    M_rot_L=zeros(3,1);
    alpha_L=0;
    U_alpha_L_dot=0;    
    for i=1:N_r-1
        r=rs(i);
        c=polyval(INSECT.cr_poly,r*1e2)*1e-2;
        xLE = polyval(INSECT.cr_LE_poly,r*1e2)*1e-2;
        xTE = polyval(INSECT.cr_TE_poly,r*1e2)*1e-2;
        xCProt = (xLE - 0.75 * c) * e1 ;
        xPivot = 0 * e1; % pitch rotation at the root
        C_rot = pi*(xLE/c - 0.75  - 0);
        
        % right wing
        U_R = r*Euler_R_dot;
        U_R_norm = U_R(1);
        omega_dot = Euler_R_dot(2);        
        sgn = 1;
        if omega_dot < 0
            sgn = -1;
        end
        
        F_rot = INSECT.rho * C_rot * U_R_norm * omega_dot * e3 * c^2 * dr * sgn  ; %* morp_factor; 
        F_rot_R = F_rot_R + F_rot ;
        M_rot_R = M_rot_R + r * cross(F_rot, e2) *-sign(U_R_norm) *sign(U_R_norm) ;
        M_rot_R = M_rot_R - cross(F_rot, 0.75*c*e1)*sign(-U_R(1))*sign(omega_dot) ;
        
        % left wing
        U_L = r*Euler_L_dot;
        U_L_norm = U_L(1); 
        alpha_L = acos(U_L_norm);
        omega_dot = Euler_L_dot(2);
        sgn = 1;
        if omega_dot < 0
            sgn = -1;
        end

        F_rot = INSECT.rho * C_rot * U_L_norm * omega_dot * e3 * c^2 * dr * sgn ; %* morp_factor;
        F_rot_L = F_rot + F_rot_L;
        M_rot_L = M_rot_L - r * cross(F_rot,e2) * -sign(U_L_norm)*sign(U_L_norm);
        M_rot_L = M_rot_L - cross(F_rot, 0.75*c*e1)*sign(-U_L(1))*sign(omega_dot) ;
     end    

    
%% CK: rotational force calculation - integrated - begin
%     hat_x_0=0.0; % CK changed to 0.25    
%     C_rot = pi*(0.75 - hat_x_0);
%     U_R = INSECT.tilde_r_2*INSECT.l*Euler_R_dot;
%     U_R_norm = U_R(1);
%     omega_dot = Euler_R_dot(2);
%     sgn = 1;
%     if U_R_norm > 0 && omega_dot > 0
%         sgn = -1;
%     elseif U_R_norm < 0 && omega_dot > 0 
%         sgn = -1;
%     end
%     F_rot_R = INSECT.rho * C_rot * U_R_norm * omega_dot * e3 * INSECT.c_bar^2 * INSECT.l * sgn  ; %* morp_factor; 
%     M_rot_R = INSECT.r_rot * cross(F_rot_R, e2) *-sign(U_R_norm) *sign(U_R_norm) ;  
%     M_rot_R = M_rot_R - cross(F_rot_R, 0.75*INSECT.c_bar*e1)*sign(-U_R(1))*sign(omega_dot) ;
%     
%     U_L = INSECT.tilde_r_2*INSECT.l*Euler_L_dot;
%     U_L_norm = U_L(1); 
%     omega_dot = Euler_L_dot(2);
%     sgn = 1;
%     if U_L_norm > 0 && omega_dot > 0
%         sgn = -1;
%     elseif U_L_norm < 0 && omega_dot > 0 
%        sgn = -1;
%     end
%     F_rot_L = INSECT.rho * C_rot * U_L_norm * omega_dot * e3 * INSECT.c_bar^2 * INSECT.l * sgn ; %* morp_factor;
%     M_rot_L = -INSECT.r_rot * cross(F_rot_L,e2) *-sign(U_L_norm)*sign(U_L_norm);
%     M_rot_L = M_rot_L;
%     M_rot_L = M_rot_L + cross(F_rot_L, 0.75*INSECT.c_bar*e1)*sign(-U_L(1))*sign(omega_dot) ;
%   
%  CK: rotational force calculation - integrated - end

end


end
%%
function [L D M alpha]=compute_LD(INSECT, U)
global e1 e2 e3

alpha=compute_alpha(U);
[C_L C_D]=wing_QS_LD_coeff(alpha);

L = 0.5 * INSECT.rho * C_L * sign(U(1)*U(3)) * cross(e2,U) * norm(U) * INSECT.S;
D = - 0.5 * INSECT.rho * C_D * norm(U) * U * INSECT.S;
M = INSECT.r_cp * cross(e2, L+D);
end

%%
function [F_rot M_rot alpha U_alpha_dot]=compute_rotational_force(INSECT, U, U_dot)
global e1 e2 e3
[alpha U_alpha_dot]=compute_alpha(U, U_dot);

%alpha_dot = compute_alpha_dot(U, U_dot);
%sgn_rot = sign(alpha)*-sign(U(3))*sign(U_alpha_dot) + (1-sign(alpha))*(-sign(U_dot(3)));
sgn_rot = sign(alpha)* -sign(U(3))*sign(U_alpha_dot) + (1-sign(alpha))*(-sign(U_dot(3)));

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight,"
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.
% theoretical value of the rotational coefficient
hat_x_0=0.25; % CK changed to 0.25

C_rot = pi*(0.75 - hat_x_0);

% From Pohly et al. Fluids 2019
%alpha_dot = compute_alpha_dot(U, U_dot);
%hat_omega = alpha_dot * c / norm(U);
%C_rot = (-11.77*hat_omega + 0.8152)*(0.75 - hat_x_0);

morp_factor = INSECT.l * INSECT.c_bar^2 * INSECT.tilde_v * INSECT.tilde_r_v_1 / INSECT.tilde_r_2;
F_rot = INSECT.rho * C_rot * abs(U_alpha_dot) * sgn_rot * e3 * morp_factor; % CK removed abs
M_rot = INSECT.r_rot * cross(e2, F_rot);

end

%%
function [F_rot M_rot alpha U_alpha_dot]=compute_rotational_force_blade(INSECT, xLE, U, U_dot, c)
global e1 e2 e3
[alpha U_alpha_dot]=compute_alpha(U, U_dot);

%alpha_dot = compute_alpha_dot(U, U_dot);
%sgn_rot = sign(alpha)*-sign(U(3))*sign(U_alpha_dot) + (1-sign(alpha))*(-sign(U_dot(3)));
sgn_rot = sign(alpha)* -sign(U(3))*sign(U_alpha_dot) + (1-sign(alpha))*(-sign(U_dot(3)));

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight,"
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.
% theoretical value of the rotational coefficient
hat_x_0= xLE/c;  %0.25; % CK changed to 0.25

C_rot = pi*(0.75 - hat_x_0);

% From Pohly et al. Fluids 2019
%alpha_dot = compute_alpha_dot(U, U_dot);
%hat_omega = alpha_dot * c / norm(U);
%C_rot = (-11.77*hat_omega + 0.8152)*(0.75 - hat_x_0);

morp_factor = INSECT.l * INSECT.c_bar^2 * INSECT.tilde_v * INSECT.tilde_r_v_1 / INSECT.tilde_r_2;
F_rot = INSECT.rho * C_rot * (U_alpha_dot) * sgn_rot * e3 ; %* morp_factor; % CK removed abs
%M_rot = INSECT.r_rot * cross(e2, F_rot);
M_rot = 0;

end

%%
function [alpha varargout]=compute_alpha(U, varargin)
%compute_alpha: the angle of attack
% [alpha]=compute_alpha(U) computes the angle of attack in [0,pi/2]
% [alpha, U_alpha_dot]=compute_alpha(U,U_dot) computes the angle of attack
% and the time-derivative of the angle of attack multiplied by the norm of
% U

norm_U = norm(U);

if norm_U == 0
    alpha=0;
    varargout{1}=0;
else
    u=U/norm_U; % unit-vector along U
    %alpha=acos(abs(u(1)));
    alpha=acos(u(1));
    %alpha=asin(u(3));
    %if acos(abs(u(1))) ~= acos(u(1))
    %    error([num2str(acos(abs(u(1)))),', ',num2str(acos(u(1)))]);
    %end
    
    if nargin > 1
        U_dot=varargin{1};
        
        if alpha > pi/4
            U_alpha_dot = 1/sin(alpha)* (cos(alpha)*dot(u,U_dot) - sign(U(1))*U_dot(1));
        else
            %U_alpha_dot = 1/sin(alpha)* (cos(alpha)*dot(u,U_dot) - sign(U(1))*U_dot(1));
            U_alpha_dot = 1/cos(alpha)* (-sin(alpha)*dot(u,U_dot)+sign(U(3))*U_dot(3));
        end
        
        varargout{1}=U_alpha_dot;
    end
end

end



%%
function [C_L C_D] = wing_QS_LD_coeff(alpha)
%wing_QS_LD_coeff: compute C_L and C_D
%[C_L C_D] = trans_force_coeff(alpha) computes the lift coefficient and the
%drag coefficient for a given angle of attack, alpha in RADIAN

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight," 
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.

% convert radian to degree
alpha_deg=alpha*180/pi;

C_L = 0.225 + 1.58 * sind(2.13*alpha_deg -7.2);
C_D = 1.92 - 1.55 * cosd(2.04*alpha_deg-9.82);
end

