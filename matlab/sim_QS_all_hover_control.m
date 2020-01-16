function sim_QS_all_hover_control
% simulate the complete dynamics (x,R,Q_R,Q_L,Q_A) for given torque acting
% on the joint to regenerate a hover trajectory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data');
des = load('sim_QS_x_hover.mat',...
    'INSECT', 't', 'N', 'x', 'x_dot', 'R', 'Q_R', 'Q_L', 'Euler_R','W_R', 'W_L', 'f_tau',...
    'x0', 'x_dot0', 'Q_A', 'W', 'W_A', 'WK', 'Euler_R_dot', 'tau');
filename='sim_QS_all_hover_control';
INSECT = des.INSECT;
WK = des.WK;

des.Euler_R_m = mean(des.Euler_R, 2);
des.Euler_L_m = mean(des.Euler_R, 2);
des.W_R_m = mean(des.W_R, 2);
des.W_L_m = mean(des.W_L, 2);
des.f_tau_fit = cell(15, 1);
for i=1:15
    des.f_tau_fit{i} = fit(des.t, des.f_tau(i, :)', 'fourier8');
end
gains.Kp_pos = -2;
gains.Kd_pos = 2;
gains.Kp_att = 0;
gains.Kd_att = 0;
%%% These two don't match; so there are inaccuracies?
% plot(t*WK.f, f_tau(5, :)); hold on; y = -(tau(4:6, :) + tau(7:9, :) + tau(10:12, :)); plot(t*WK.f, y(2, :), 'r');

eps = 0;
x0 = des.x0 + rand(3,1)*eps;
R0 = des.R(:, :, 1)*expmso3(rand(3,1)*eps);
Q_R0 = des.Q_R(:, :, 1)*expmso3(rand(3,1)*eps);
Q_L0 = des.Q_L(:, :, 1)*expmso3(rand(3,1)*eps);
Q_A0 = des.Q_A(:, :, 1)*expmso3(rand(3,1)*eps);
x_dot0 = des.x_dot0 + rand(3,1)*eps;
W0 = des.W(:, 1) + rand(3,1)*eps;
W_R0 = des.W_R(:, 1) + rand(3,1)*eps;
W_L0 = des.W_L(:, 1) + rand(3,1)*eps;
W_A0 = des.W_A(:, 1) + rand(3,1)*eps;

X0=[x0; reshape(R0,9,1); reshape(Q_R0,9,1); reshape(Q_L0,9,1); reshape(Q_A0,9,1);...
    x_dot0; W0; W_R0; W_L0; W_A0];

N = 1001;
t = linspace(0, 2/WK.f, N);
dt = t(2) - t(1);
i = 1;
[t X]=ode15s(@(t,X) eom(INSECT, t, X, des, gains, i, dt), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

% N = 5001;
% t = linspace(0, 2/WK.f, N);
% X = zeros(N, 54);
% X(1, :) = X0;
% dt = t(2) - t(1);
% % Implicit trapezoidal method
% options = optimoptions('fsolve', 'Display', 'none', 'UseParallel', true);
% for i=1:(N-1)
%     X_i = X(i, :)';
%     f_i = eom(INSECT, t(i), X_i, des, gains, i, dt);
%     fun = @(X) (X - X_i - 0.5 * dt *(f_i +...
%         eom(INSECT, t(i) + dt, X, des, gains, i, dt)));
%     X_init = X_i + dt * f_i;
%     X(i+1, :) = fsolve(fun, X_init, options)';
% end

x=X(:,1:3)';
x_dot=X(:,40:42)';
W=X(:,43:45)';
W_R=X(:,46:48)';
W_L=X(:,49:51)';
W_A=X(:,52:54)';

R=zeros(3,3,N);Q_R=zeros(3,3,N);Q_L=zeros(3,3,N);Q_A=zeros(3,3,N);
for k=1:N
    R(:,:,k)=reshape(X(k,4:12),3,3);
    Q_R(:,:,k)=reshape(X(k,13:21),3,3);
    Q_L(:,:,k)=reshape(X(k,22:30),3,3);
    Q_A(:,:,k)=reshape(X(k,31:39),3,3);
end

save(filename);
evalin('base',['load ' filename]);
end

function [X_dot F_R F_L] = eom(INSECT, t, X, des, gains, i, dt)
x=X(1:3);
R=reshape(X(4:12),3,3);
Q_R=reshape(X(13:21),3,3);
Q_L=reshape(X(22:30),3,3);
Q_A=reshape(X(31:39),3,3);
x_dot=X(40:42);
W=X(43:45);
W_R=X(46:48);
W_L=X(49:51);
W_A=X(52:54);

xi=[x_dot; W; W_R; W_L; W_A];
[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

[L_R L_L D_R D_L M_R M_L ...
    F_rot_R F_rot_L M_rot_R M_rot_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, zeros(3, 1), zeros(3, 1), x_dot, R, W, Q_R, Q_L);
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
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);

% [L_R_des, L_L_des]=wing_QS_aerodynamics(INSECT, des.W_R_m, des.W_L_m, zeros(3, 1), zeros(3, 1));
% K_R = -L_R_des(1)/des.W_R_m(1);
% K_L = -L_L_des(1)/des.W_L_m(1);
% 
% d_x = des.x(:, i) - x;
% d_x_dot = des.x_dot(:, i) - x_dot;
% pos_err = INSECT.m*(gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot) /2;
% d_att_R = (1/K_R * des.Q_R_m' * des.R_m' * pos_err);
% Euler_R = des.Euler_R_m;
% T_R = [-cos(Euler_R(3))*cos(Euler_R(2)), -sin(Euler_R(2));
%        -cos(Euler_R(3))*sin(Euler_R(2)), cos(Euler_R(2))];
% d_att_R = T_R \ [d_att_R(1); d_att_R(3);];
% d_att_L = (1/K_L * des.Q_L_m' * des.R_m' * pos_err);
% Euler_L = des.Euler_L_m;
% T_L = [cos(Euler_L(3))*cos(Euler_L(2)), sin(Euler_L(2));
%        cos(Euler_L(3))*sin(Euler_L(2)), -cos(Euler_L(2))];
% d_att_L = T_L \ [d_att_L(1); d_att_L(3);];
% 
% d_omega_R = des.W_R(:, i) - W_R;
% d_tau_R = des.Q_R_m * INSECT.J_R *...
%     (gains.Kp_att * [d_att_R(1); 0; d_att_R(2);] + gains.Kd_att * d_omega_R);
% d_omega_L = des.W_L(:, i) - W_L;
% d_tau_L = des.Q_L_m * INSECT.J_L *...
%     (gains.Kp_att * [d_att_L(1); 0; d_att_L(2);] + gains.Kd_att * d_omega_L);
% 
% f_tau = des.f_tau(:, i) + [zeros(6, 1); d_tau_R; d_tau_L; zeros(3, 1)];

f_tau = zeros(15, 1);
for j=4:15
    f_tau(j) = des.f_tau_fit{j}(t);
end
xi_dot=JJ\(-LL*xi + co_ad*JJ*xi - dU + f_a + f_tau);

ke = 0;
I = eye(3);
R_dot = R*hat(W) - ke*R*(R'*R - I);
Q_R_dot = Q_R*hat(W_R) - ke*Q_R*(Q_R'*Q_R - I);
Q_L_dot = Q_L*hat(W_L) - ke*Q_L*(Q_L'*Q_L - I);
Q_A_dot = Q_A*hat(W_A) - ke*Q_A*(Q_A'*Q_A - I);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_R_dot,9,1); reshape(Q_L_dot,9,1); reshape(Q_A_dot,9,1); xi_dot];
end
