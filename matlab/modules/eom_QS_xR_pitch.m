function [X_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau, Euler_R, Euler_R_dot, f_t_2]= eom_QS_xR_pitch(INSECT, WK_R, WK_L, t, X)
% Returns the states and forces for given wing kinematics and
% abdomen attitude.

x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);
theta_R=X(19);
theta_L=X(20);
theta_R_dot=X(21);
theta_L_dot=X(22);

% wing/abdomen attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
Euler_R(2) = theta_R; Euler_R_dot(2) = theta_R_dot; Euler_R_ddot(2) = 0;
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
Euler_L(2) = theta_L; Euler_L_dot(2) = theta_L_dot; Euler_L_ddot(2) = 0;
[Q_R, Q_L, W_R, W_L] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot);
[Q_A, W_A, W_A_dot, theta_A] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen

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
% T_R_second = T_R_theta*Euler_R_dot*Euler_R_dot(2) + T_R_psi*Euler_R_dot*Euler_R_dot(3);

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
% T_L_second = T_L_theta*Euler_L_dot*Euler_L_dot(2) + T_L_psi*Euler_L_dot*Euler_L_dot(3);

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

%%
% 2nd component of Euler_R_ddot has been initialized to be 0
xi_2_dot_pp = [T_R*Euler_R_ddot + (T_R_theta*Euler_R_dot(2) + T_R_psi*Euler_R_dot(3))*Euler_R_dot;
               T_L*Euler_L_ddot + (T_L_theta*Euler_L_dot(2) + T_L_psi*Euler_L_dot(3))*Euler_L_dot;
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
xi_1_dot = soln(1:6);
Euler_R_ddot(2) = soln(7);
Euler_L_ddot(2) = soln(8);

%%
W_R_dot = T_R*Euler_R_ddot + (T_R_theta*Euler_R_dot(2) + T_R_psi*Euler_R_dot(3))*Euler_R_dot;
W_L_dot = T_L*Euler_L_ddot + (T_L_theta*Euler_L_dot(2) + T_L_psi*Euler_L_dot(3))*Euler_L_dot;
xi_2_dot = [W_R_dot; W_L_dot; W_A_dot];

f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1+JJ_22*xi_2) ...
    + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2 - f_t_2;
f_tau_1 = C*f_tau_2;
f_tau = [f_tau_1; f_tau_2];

% tau = blkdiag(Q_R, Q_L, Q_A)*f_tau_2;
tau = zeros(9, 9);
tau(1:3, 1:3) = Q_R; tau(4:6, 4:6) = Q_L; tau(7:9, 7:9) = Q_A;
tau = tau*f_tau_2;

R_dot = R*hat(W);
X_dot=[x_dot; reshape(R_dot,9,1); xi_1_dot; ...
    theta_R_dot; theta_L_dot; Euler_R_ddot(2); Euler_L_ddot(2)];

end
