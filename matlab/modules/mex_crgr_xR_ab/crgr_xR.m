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
s = length(b);
x_dot = zeros(3, s);
W = zeros(3, s);
K_xi = zeros(6, s);

x1 = X(1:3);
R1 = reshape(X(4:12),3,3);
x_dot(:, 1) = X(13:15);
W(:, 1) = X(16:18);

xi1 = X(13:18);
t1 = t+c(1)*dt;
K_xi(:, 1) = deriv_xi(INSECT, x1, R1, x_dot(:, 1), W(:, 1), WK_R, WK_L, t1);

for i = 2:s
    x = x1 + dt* sum(bsxfun(@times, a(i,1:(i-1)), x_dot(:, 1:(i-1))), 2);
    R = R1;
    for j = 1:(i-1)
        R = R * expmhat(dt*a(i, j)*W(:, j));
    end
    xi = xi1 + dt* sum(bsxfun(@times, a(i,1:(i-1)), K_xi(:, 1:(i-1))), 2);
    x_dot(:, i) = xi(1:3); W(:, i) = xi(4:6);
    ti = t+c(i)*dt;
    K_xi(:, i) = deriv_xi(INSECT, x, R, x_dot(:, i), W(:, i), WK_R, WK_L, ti);
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
