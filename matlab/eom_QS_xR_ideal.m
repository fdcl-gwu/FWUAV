function [X_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau, Euler_R, Euler_R_dot]= eom_QS_xR_ideal(INSECT, WK_R, WK_L, t, X, u_control, tmp_f_des)
% Returns the states and forces for given wing kinematics and
% abdomen attitude.

x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);

% wing/abdomen attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A, W_A, W_A_dot, theta_A] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen

[L_R, L_L, D_R, D_L, M_R, M_L, ...
    F_rot_R, F_rot_L, M_rot_R, M_rot_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
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
%f_g=zeros(15,1);

% Euler-Lagrange equation
xi_1=[x_dot; W];
xi_2=[W_R; W_L; W_A];
xi_2_dot=[W_R_dot; W_L_dot; W_A_dot];

[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
% co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));
co_ad=zeros(15, 15);
co_ad(4:6, 4:6) = -hat(W); co_ad(7:9, 7:9) = -hat(W_R);
co_ad(10:12, 10:12) = -hat(W_L); co_ad(13:15, 13:15) = -hat(W_A);

[JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_6_9(JJ);
[LL_11, LL_12, LL_21, LL_22] = inertia_sub_decompose_6_9(LL);
[co_ad_11, ~, ~, co_ad_22] = inertia_sub_decompose_6_9(co_ad);

C=[zeros(3,9);
    -Q_R -Q_L -Q_A];

tmp_1 = -(co_ad_11*JJ_11-C*co_ad_22*JJ_21)*xi_1 + (LL_11-C*LL_21)*xi_1;
tmp_2 = -(JJ_12-C*JJ_22)*xi_2_dot + (co_ad_11*JJ_12-C*co_ad_22*JJ_22)*xi_2 ...
    -(LL_12-C*LL_22)*xi_2;
tmp_f = f_a_1+f_g_1 - C*(f_a_2+f_g_2);

% xi_1_dot=(JJ_11-C*JJ_21)\(-tmp_1+tmp_2+tmp_f+u_control);
xi_1_dot=(JJ_11-C*JJ_21)\(-tmp_1+tmp_2+tmp_f_des+u_control);
% xi_1_dot=(JJ_11-C*JJ_21)\(tmp_2+tmp_f_des+u_control);
% xi_1_dot=(JJ_11-C*JJ_21)\(u_control);

f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1+JJ_22*xi_2) ...
    + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
f_tau_1 = C*f_tau_2;
f_tau = [f_tau_1; f_tau_2];

% tau = blkdiag(Q_R, Q_L, Q_A)*f_tau_2;
tau = zeros(9, 9);
tau(1:3, 1:3) = Q_R; tau(4:6, 4:6) = Q_L; tau(7:9, 7:9) = Q_A;
tau = tau*f_tau_2;

% xi=[xi_1;xi_2];
% xi_dot=JJ\( co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau);
% disp(norm(xi_dot - [xi_1_dot; xi_2_dot]));

R_dot = R*hat(W);
X_dot=[x_dot; reshape(R_dot,9,1); xi_1_dot];

end
