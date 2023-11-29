function X = rk_xR_control(INSECT, WK_R_des, WK_L_des, t, X0, dang0, param, N_dang, N_con, N_per_iter)
%% assuming equally spaced t array
N = length(t);
dt = t(2) - t(1);
X = zeros(N, length(X0));
X(1, :) = X0;

% Jackiewicz, Z., Marthinsen, A. and Owren, B., 2000.
% Construction of Runge–Kutta methods of Crouch–Grossman type of high order.
% Advances in Computational Mathematics, 13(4), pp.405-415.

% 4th order Kutta
a = [0, 0, 0, 0;
    1/3, 0, 0, 0;
    -1/3, 1, 0, 0;
    1, -1, 1, 0];
b = [1/8, 3/8, 3/8, 1/8];
c = [0, 1/3, 2/3, 1];

% ode45 Dormand–Prince method
% a = , b = , c = 

for con=1:N_con
    param_idx = (1+(con-1)*N_dang):(con*N_dang);
    param_con = param(param_idx);
    idx_1 = 1+(con-1)*N_per_iter; idx_end = 1+con*N_per_iter;
    
    t0 = t(idx_1);
    for k=idx_1:(idx_end-1)
        X(k+1, :) = runge_kutta(INSECT, X(k, :)', t(k), dt, a, b, c, WK_R_des, WK_L_des, t0, dang0, param_con);
    end
    dang0 = dang0 + param_con * (t(idx_end) - t0);
end

end

function Xout = runge_kutta(INSECT, X, t, dt, a, b, c, WK_R_des, WK_L_des, t0, dang0, param_con)
%% 4th order Crouch Grossman
s = length(b);
K_xi = zeros(18, s);
xi1 = X(1:18);

for i = 1:s
    xi = xi1 + dt* sum(bsxfun(@times, a(i,1:(i-1)), K_xi(:, 1:(i-1))), 2);
    x = xi(1:3); R = reshape(xi(4:12), 3, 3);
    x_dot = xi(13:15); W = xi(16:18);
    ti = t+c(i)*dt;
    dang = dang0 + param_con*(ti - t0);
    [WK_R, WK_L] = get_WK(WK_R_des, WK_L_des, dang);
    K_xi(:, i) = deriv_full(INSECT, x, R, x_dot, W, WK_R, WK_L, ti);
end

% Final values
xiout = xi1 + dt* sum(bsxfun(@times, b, K_xi), 2);
% xiout = xi1 + dt* K_xi * b'; % little slower

Xout = xiout;

end

function [WKR_new, WKL_new] = get_WK(WKR_old, WKL_old, dang)
%%
WKR_new = WKR_old; WKL_new = WKL_old;
WKR_new.phi_m = WKR_new.phi_m + dang(1) + dang(3);
WKL_new.phi_m = WKL_new.phi_m + dang(1) - dang(3);
WKR_new.theta_0 = WKR_new.theta_0 + dang(2) + dang(5);
WKL_new.theta_0 = WKL_new.theta_0 + dang(2) - dang(5);
WKR_new.phi_0 = WKR_new.phi_0 + dang(4);
WKL_new.phi_0 = WKL_new.phi_0 + dang(4);
WKR_new.psi_m = WKR_new.psi_m + dang(6);
WKL_new.psi_m = WKL_new.psi_m - dang(6);
    
% Include more wing parameters?

%     WKR_new.psi_m = WKR_new.psi_m + dang(4) + dang(6);
%     WKL_new.psi_m = WKL_new.psi_m + dang(4) - dang(6);
    
%     WKR_new.theta_A_m = WKR_new.theta_A_m + dang(7);
%     WKL_new.theta_A_m = WKL_new.theta_A_m + dang(7);
end

function xi_dot = deriv_full(INSECT, x, R, x_dot, W, WK_R, WK_L, t)
    xi_dot = [x_dot;
              reshape(R*hat(W), 9, 1);
              deriv_xi(INSECT, x, R, x_dot, W, WK_R, WK_L, t)];
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

