function [L_R L_L D_R D_L M_R M_L ...
    F_rot_R F_rot_L M_rot_R M_rot_L ...
    alpha_R alpha_L U_alpha_R_dot U_alpha_L_dot U_R U_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, varargin)
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

if nargin <6
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
    
else
    %% integrate infinitesimal forces over the wing span
    x_dot=varargin{1};
    R=varargin{2};
    W=varargin{3};
    Q_R=varargin{4};
    Q_L=varargin{5};
    
    N_r=50;
    rs=linspace(0,INSECT.l,N_r);
    dr=rs(2)-rs(1);
    L_R = zeros(3,1);
    D_R = zeros(3,1);
    M_R = zeros(3,1);
    L_L = zeros(3,1);
    D_L = zeros(3,1);
    M_L = zeros(3,1);
    
    for i=1:N_r-1
        r=rs(i);
        c=polyval(INSECT.cr_poly,r/INSECT.scale)*INSECT.scale;
        U_R = (eye(3)-e2*e2')*Q_R'*(R'*x_dot+hat(W)*INSECT.mu_R) + r*hat(Q_R*W+W_R)*e2;
        [L D M]=compute_LD(INSECT, U_R);
        L_R = L_R + L/INSECT.S*c*dr;
        D_R = D_R + D/INSECT.S*c*dr;
        M_R = M_R + cross(r*e2,(L+D))/INSECT.S*c*dr;
        
        U_L = (eye(3)-e2*e2')*Q_L'*(R'*x_dot+hat(W)*INSECT.mu_L) - r*hat(Q_L*W+W_L)*e2;
        [L D M]=compute_LD(INSECT, U_L);
        L_L = L_L + L/INSECT.S*c*dr;
        D_L = D_L + D/INSECT.S*c*dr;
        M_L = M_L + cross(-r*e2,(L+D))/INSECT.S*c*dr;
    end
    
    % rotational force and moment are set to zero
    F_rot_R=zeros(3,1);
    M_rot_R=zeros(3,1);
    alpha_R=0;
    U_alpha_R_dot=0;
    
    F_rot_L=zeros(3,1);
    M_rot_L=zeros(3,1);
    alpha_L=0;
    U_alpha_L_dot=0;
end


end

function [L D M alpha]=compute_LD(INSECT, U)
global e1 e2 e3

alpha=compute_alpha(U);
[C_L C_D]=wing_QS_LD_coeff(INSECT, alpha);

L = 0.5 * INSECT.rho * C_L * sign(U(1)*U(3)) * cross(e2,U) * norm(U) * INSECT.S;
D = - 0.5 * INSECT.rho * C_D * norm(U) * U * INSECT.S;
M = INSECT.r_cp * cross(e2, L+D);
end

function [F_rot M_rot alpha U_alpha_dot]=compute_rotational_force(INSECT, U, U_dot)
global e1 e2 e3
[alpha U_alpha_dot]=compute_alpha(U, U_dot);

sgn_rot = sign(alpha)*-sign(U(3))*sign(U_alpha_dot) + (1-sign(alpha))*(-sign(U_dot(3)));

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight,"
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.
% theoretical value of the rotational coefficient
hat_x_0=0.5;
C_rot = pi*(0.75 - hat_x_0);
morp_factor = INSECT.l * INSECT.c_bar^2 * INSECT.tilde_v * INSECT.tilde_r_v_1 / INSECT.tilde_r_2;
F_rot = INSECT.rho * C_rot * abs(U_alpha_dot) * sgn_rot * e3 * morp_factor;
M_rot = INSECT.r_rot * cross(e2, F_rot);

end

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
    alpha=acos(abs(u(1)));
    
    if nargin > 1
        U_dot=varargin{1};
        
        if alpha > pi/4
            U_alpha_dot = 1/sin(alpha)* (cos(alpha)*dot(u,U_dot)-sign(U(1))*U_dot(1));
        else
            U_alpha_dot = 1/cos(alpha)* (-sin(alpha)*dot(u,U_dot)+sign(U(3))*U_dot(3));
        end
        
        varargout{1}=U_alpha_dot;
    end
end

end

function [C_L C_D] = wing_QS_LD_coeff(INSECT, alpha)
%wing_QS_LD_coeff: compute C_L and C_D
%[C_L C_D] = trans_force_coeff(alpha) computes the lift coefficient and the
%drag coefficient for a given angle of attack, alpha in RADIAN

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight," 
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.

% Berman, Gordon J., and Z. Jane Wang. "Energy-minimizing kinematics in
% hovering insect flight." Journal of Fluid Mechanics 582 (2007): 153-168.

% convert radian to degree
alpha_deg=alpha*180/pi;

switch INSECT.name
    case 'MONARCH'
        C_L = 0.225 + 1.58 * sind(2.13*alpha_deg -7.2);
        C_D = 1.92 - 1.55 * cosd(2.04*alpha_deg-9.82);
    otherwise
        C_L = INSECT.C_T * sin(2*alpha);
        C_D = INSECT.C_D_0 * cos(alpha)^2 + INSECT.C_D_pi2 * sin(alpha)^2;
end

end
