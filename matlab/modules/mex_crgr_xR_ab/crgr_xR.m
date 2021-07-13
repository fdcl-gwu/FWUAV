function X = crgr_xR(INSECT, WK_R, WK_L, t, X0)
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
x1=X(1:3);
R1=reshape(X(4:12),3,3);
x_dot1=X(13:15);
W1=X(16:18);

xi1=X(13:18);
t1=t+c(1)*dt;
K1_xi = deriv_xi(INSECT, x1, R1, x_dot1, W1, WK_R, WK_L, t1);

x2 = x1 + dt*a(2,1)* x_dot1;
R2 = R1 * expmhat(dt*a(2,1)*W1);
xi2 = xi1 + dt*a(2,1)* K1_xi;
x_dot2 = xi2(1:3); W2 = xi2(4:6);
t2=t+c(2)*dt;
K2_xi = deriv_xi(INSECT, x2, R2, x_dot2, W2, WK_R, WK_L, t2);

x3 = x1 + dt* (a(3,1)* x_dot1 + a(3,2)* x_dot2);
R3 = R1 * expmhat(dt*a(3,1)*W1) * expmhat(dt*a(3,2)*W2);
xi3 = xi1 + dt* (a(3,1)* K1_xi + a(3,2)* K2_xi);
x_dot3 = xi3(1:3); W3 = xi3(4:6);
t3=t+c(3)*dt;
K3_xi = deriv_xi(INSECT, x3, R3, x_dot3, W3, WK_R, WK_L, t3);

x4 = x1 + dt* (a(4,1)* x_dot1 + a(4,2)* x_dot2 + a(4,3)* x_dot3);
R4 = R1 * expmhat(dt*a(4,1)*W1) * expmhat(dt*a(4,2)*W2) * expmhat(dt*a(4,3)*W3);
xi4 = xi1 + dt* (a(4,1)* K1_xi + a(4,2)* K2_xi + a(4,3)* K3_xi);
x_dot4 = xi4(1:3); W4 = xi4(4:6);
t4=t+c(4)*dt;
K4_xi = deriv_xi(INSECT, x4, R4, x_dot4, W4, WK_R, WK_L, t4);

x5 = x1 + dt* (a(5,1)* x_dot1 + a(5,2)* x_dot2 + a(5,3)* x_dot3 + a(5,4)* x_dot4);
R5 = R1 * expmhat(dt*a(5,1)*W1) * expmhat(dt*a(5,2)*W2) * expmhat(dt*a(5,3)*W3) * expmhat(dt*a(5,4)*W4);
xi5 = xi1 + dt* (a(5,1)* K1_xi + a(5,2)* K2_xi + a(5,3)* K3_xi + a(5,4)* K4_xi);
x_dot5 = xi5(1:3); W5 = xi5(4:6);
t5=t+c(5)*dt;
K5_xi = deriv_xi(INSECT, x5, R5, x_dot5, W5, WK_R, WK_L, t5);

xout = x1 + dt* (b(1)* x_dot1 + b(2)* x_dot2 + b(3)* x_dot3 + b(4)* x_dot4 + b(5)* x_dot5);
Rout = R1 * expmhat(dt*b(1)*W1) * expmhat(dt*b(2)*W2) * expmhat(dt*b(3)*W3) * expmhat(dt*b(4)*W4) * expmhat(dt*b(5)*W5);
xiout = xi1 + dt* (b(1)* K1_xi + b(2)* K2_xi + b(3)* K3_xi + b(4)* K4_xi + b(5)* K5_xi);

Xout = [xout; reshape(Rout, 9, 1); xiout];

end

function xi_1_dot = deriv_xi(INSECT, x, R, x_dot, W, WK_R, WK_L, t)
%%
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen

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

xi_1_dot=(JJ_11-C*JJ_21)\(-tmp_1+tmp_2+tmp_f);

end
