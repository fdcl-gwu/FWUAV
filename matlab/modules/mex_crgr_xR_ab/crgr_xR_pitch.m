function X = crgr_xR_pitch(INSECT, WK_R, WK_L, t, X0)
%% assuming equally spaced t array
N = length(t);
dt = t(2) - t(1);
X = zeros(N, length(X0));
X(1, :) = X0;

% Jackiewicz, Z., Marthinsen, A. and Owren, B., 2000.
% Construction of Runge–Kutta methods of Crouch–Grossman type of high order.
% Advances in Computational Mathematics, 13(4), pp.405-415.
a = [0, 0, 0, 0, 0;
    0.8177227988124852, 0, 0, 0, 0;
    0.3199876375476427, 0.0659864263556022, 0, 0, 0;
    0.9214417194464946, 0.4997857776773573, -1.0969984448371582, 0, 0;
    0.3552358559023322, 0.2390958372307326, 1.3918565724203246, -1.1092979392113465, 0];
b = [0.1370831520630755, -0.0183698531564020, 0.7397813985370780, -0.1907142565505889, 0.3322195591068374];
c = [0, 0.8177227988124852, 0.3859740639032449, 0.3242290522866937, 0.8768903263420429];

for k=1:N-1
    X(k+1, :) = crouch_grossman_4th(INSECT, X(k, :)', t(k), dt, a, b, c, WK_R, WK_L);
end

end

function Xout = crouch_grossman_4th(INSECT, X, t, dt, a, b, c, WK_R, WK_L)
%% 4th order Crouch Grossman
s = length(b);
x_dot = zeros(3, s);
W = zeros(3, s);
K_xi = zeros(10, s);

x1 = X(1:3);
R1 = reshape(X(4:12),3,3);
x_dot(:, 1) = X(13:15);
W(:, 1) = X(16:18);
thetas1 = X(19:22);
xi1 = X(13:22);
t1 = t+c(1)*dt;
K_xi(:, 1) = deriv_xi(INSECT, x1, R1, x_dot(:, 1), W(:, 1), thetas1, WK_R, WK_L, t1);

for i = 2:s
    x = x1 + dt* sum(bsxfun(@times, a(i,1:(i-1)), x_dot(:, 1:(i-1))), 2);
    R = R1;
    for j = 1:(i-1)
        R = R * expmhat(dt*a(i, j)*W(:, j));
    end
    xi = xi1 + dt* sum(bsxfun(@times, a(i,1:(i-1)), K_xi(:, 1:(i-1))), 2);
    x_dot(:, i) = xi(1:3); W(:, i) = xi(4:6); thetas = xi(7:10);
    ti = t+c(i)*dt;
    K_xi(:, i) = deriv_xi(INSECT, x, R, x_dot(:, i), W(:, i), thetas, WK_R, WK_L, ti);
end

% Final values
xout = x1 + dt* sum(bsxfun(@times, b, x_dot), 2);
Rout = R1;
for j = 1:s
    Rout = Rout * expmhat(dt*b(j)*W(:, j));
end
xiout = xi1 + dt* sum(bsxfun(@times, b, K_xi), 2);
% xiout = xi1 + dt* K_xi * b'; % little slower

Xout = [xout; reshape(Rout, 9, 1); xiout];

end

function xi_1_dot = deriv_xi(INSECT, x, R, x_dot, W, thetas, WK_R, WK_L, t)
%%
theta_R = thetas(1); theta_L = thetas(2);
theta_R_dot = thetas(3); theta_L_dot = thetas(4);

% wing/abdomen attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
Euler_R(2) = theta_R; Euler_R_dot(2) = theta_R_dot; Euler_R_ddot(2) = 0;
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
Euler_L(2) = theta_L; Euler_L_dot(2) = theta_L_dot; Euler_L_ddot(2) = 0;
[Q_R, Q_L, W_R, W_L] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot);
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen

[L_R, L_L, D_R, D_L, M_R, M_L, ...
    F_rot_R, F_rot_L, M_rot_R, M_rot_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, 0, 0, x_dot, R, W, Q_R, Q_L);
F_R=L_R+D_R+F_rot_R;
F_L=L_L+D_L+F_rot_L;
M_R=M_R+M_rot_R;
M_L=M_L+M_rot_L;
M_A=zeros(3,1);

f_a=[R*Q_R*F_R + R*Q_L*F_L;
    hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R;
    M_L;
    M_A];
f_a_1=f_a(1:6);
f_a_2=f_a(7:15);
% f_a=zeros(15,1);

% gravitational force and moment
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:6);
f_g_2=f_g(7:15);

% Torsional spring at wing root
e2 = [0, 1, 0]';
f_t_2=[-(INSECT.Ktorsion*theta_R + INSECT.Ctorsion*theta_R_dot - INSECT.ftau0)*e2;
       -(INSECT.Ktorsion*theta_L + INSECT.Ctorsion*theta_L_dot - INSECT.ftau0)*e2;
        zeros(3,1)];

psi_R=Euler_R(3);
T_R=[cos(psi_R)*cos(theta_R), 0, sin(theta_R);
    sin(psi_R), 1, 0;
    cos(psi_R)*sin(theta_R), 0, -cos(theta_R)];
T_R_theta=[ -cos(psi_R)*sin(theta_R), 0, cos(theta_R);
    0, 0, 0;
    cos(psi_R)*cos(theta_R), 0, sin(theta_R)];
T_R_psi=[ -cos(theta_R)*sin(psi_R), 0, 0;
    cos(psi_R), 0, 0;
    -sin(psi_R)*sin(theta_R), 0, 0];
T_R_second = (T_R_theta*Euler_R_dot(2) + T_R_psi*Euler_R_dot(3)) * Euler_R_dot;

psi_L=Euler_L(3);
T_L = [-cos(psi_L)*cos(theta_L), 0, -sin(theta_L);
    sin(psi_L), 1, 0;
    -cos(psi_L)*sin(theta_L), 0, cos(theta_L)];
T_L_theta = [cos(psi_L)*sin(theta_L), 0, -cos(theta_L);
    0, 0, 0;
    -cos(psi_L)*cos(theta_L), 0, -sin(theta_L)];
T_L_psi = [sin(psi_L)*cos(theta_L), 0, 0;
    cos(psi_L), 0, 0;
    sin(psi_L)*sin(theta_L), 0, 0];
T_L_second = (T_L_theta*Euler_L_dot(2) + T_L_psi*Euler_L_dot(3)) * Euler_L_dot;

% Euler-Lagrange equation
xi_1=[x_dot; W];
xi_2=[W_R; W_L; W_A];

[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
% co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));
co_ad=zeros(15, 15);
co_ad(4:6, 4:6) = -hat(W); co_ad(7:9, 7:9) = -hat(W_R);
co_ad(10:12, 10:12) = -hat(W_L); co_ad(13:15, 13:15) = -hat(W_A);

[JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_6_9(JJ);
[LL_11, LL_12, LL_21, LL_22] = inertia_sub_decompose_6_9(LL);
[co_ad_11, ~, ~, co_ad_22] = inertia_sub_decompose_6_9(co_ad);

% 2nd component of Euler_R_ddot has been initialized to be 0
xi_2_dot_pp = [T_R*Euler_R_ddot + T_R_second;
               T_L*Euler_L_ddot + T_L_second;
               W_A_dot];
T_ddot = zeros(9, 2);
T_ddot(1:3, 1) = T_R(:, 2); T_ddot(4:6, 2) = T_L(:, 2);
C=[zeros(3,9);
    -Q_R -Q_L -Q_A];
tmp_1 = -(co_ad_11*JJ_11-C*co_ad_22*JJ_21)*xi_1 + (LL_11-C*LL_21)*xi_1;
tmp_2 = (co_ad_11*JJ_12-C*co_ad_22*JJ_22)*xi_2 - (LL_12-C*LL_22)*xi_2;
tmp_f = f_a_1+f_g_1 - C*(f_a_2+f_g_2+f_t_2);

eqn_matrix = zeros(8, 8); rhs = zeros(8, 1);

E_RL = zeros(2, 9); E_RL(1, 2) = 1; E_RL(2, 5) =1;
eqn_matrix(1:6, :) = [(JJ_11-C*JJ_21), (JJ_12-C*JJ_22)*T_ddot];
eqn_matrix(7:8, :) = E_RL * [JJ_21, JJ_22*T_ddot];

rhs(1:6) = -tmp_1+tmp_2+tmp_f - (JJ_12-C*JJ_22)*xi_2_dot_pp;
rhs(7:8) = E_RL * (co_ad_22*(JJ_21*xi_1 + JJ_22*xi_2) - LL_21*xi_1 - LL_22*xi_2 + ...
                   f_a_2 + f_g_2 + f_t_2 - JJ_22*xi_2_dot_pp);

soln = eqn_matrix \ rhs;
xi_1_dot = [soln(1:6); theta_R_dot; theta_L_dot; soln(7:8)];

end
