function sim_QS_x_hover_control
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
% close all;
addpath('./modules', './sim_data');
des=load('sim_QS_x_hover.mat',...
    'INSECT', 't', 'N', 'x', 'x_dot', 'R', 'Q_R', 'Q_L', 'W_R', 'W_L', 'f_tau',...
    'x0', 'x_dot0', 'Q_A', 'WK');

filename='sim_QS_x_hover_control';
INSECT = des.INSECT;
WK = des.WK;
des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1);
% fttype = "a0 + ";
% for i=1:20
%     fttype = fttype + 'a' + string(i) + '*cos(' + string(i) + '*x*w) + b'...
%         + string(i) + '*sin(' + string(i) + '*x*w) +';
% end
% fttype = fttype + "0";
% fttype = fittype(fttype);
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', 'fourier8');
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', 'fourier8');
end

% Values obtained from parametric study
% Orig values
des.df_a_1_by_dphi_m = 1.5e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0
des.df_a_1_by_dtheta_m = -1.4e-3 / 0.1; % dtheta_m_R > 0, dtheta_m_L > 0
des.df_a_2_by_dpsi_m = 1e-3 / 0.1; % dpsi_m_R > 0, dpsi_m_L < 0
des.df_a_3_by_dphi_m = 1.3e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0
% New values
des.df_a_1_by_dphi_m = 0.54e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0
des.df_a_1_by_dtheta_m = -0.6e-3 / 0.1; % dtheta_m_R > 0, dtheta_m_L > 0
des.df_a_2_by_dpsi_m = 0.3e-3 / 0.1; % dpsi_m_R > 0, dpsi_m_L < 0
des.df_a_3_by_dphi_m = 0.47e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0

% Gains
% Original
gains.Kp_pos = -2;
gains.Kd_pos = 2;
% New
gains.Kp_pos = 2;
gains.Kd_pos = 2;
gains.Ki_pos = 0;

%% Simulation
eps1 = 1e-2; eps2 = 1e-1;
x0 = des.x0 + rand(3,1)*eps1;
x_dot0 = des.x_dot0 + rand(3,1)*eps2;
X0 = [x0; x_dot0;];

N = 2001;
N_period = 20;
N_single = round((N-1)/N_period);
T = N_period/WK.f;
t = linspace(0,T,N);

% [t, X]=ode45(@(t,X) eom(INSECT, WK, WK, t, X, des, gains, i), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

X = zeros(N, 9);
X(1, :) = [X0; zeros(3, 1)];
dt = t(2) - t(1);
% % Explicit RK4
for i=1:(N-1)
%     if mod(i-1, N_single) == 0
%         x0 = X(i, 1:3)';
%     end
    k1 = dt * eom(INSECT, WK, WK, t(i), X(i, :)', des, gains, i, x0);
    k2 = dt * eom(INSECT, WK, WK, t(i)+dt/2, X(i, :)'+k1/2, des, gains, i, x0);
    k3 = dt * eom(INSECT, WK, WK, t(i)+dt/2, X(i, :)'+k2/2, des, gains, i, x0);
    k4 = dt * eom(INSECT, WK, WK, t(i), X(i, :)'+k3, des, gains, i, x0);
    X(i+1, :) = X(i, :) + 1/6 * (k1 + 2*k2 + 2*k3 + k4)';
end

x=X(:,1:3)';
x_dot=X(:,4:6)';

R=zeros(3,3,N);
for k=1:N
%     if mod(k-1, N_single) == 0
%         x0 = X(k, 1:3)';
%     end
    [X_dot(:,k), R(:,:,k) Q_R(:,:,k) Q_L(:,:,k) Q_A(:,:,k) theta_B(k) theta_A(k) W(:,k) W_dot(:,k) W_R(:,k) W_R_dot(:,k) W_L(:,k) W_L_dot(:,k) W_A(:,k) W_A_dot(:,k) F_R(:,k) F_L(:,k) M_R(:,k) M_L(:,k) f_a(:,k) f_g(:,k) f_tau(:,k) tau(:,k) Euler_R(:,k) Euler_R_dot(:,k) pos_err(:, k)]= eom(INSECT, WK, WK, t(k), X(k,:)', des, gains, k, x0);
    F_B(:,k)=Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);
end
x_ddot = X_dot(4:6,:);

%% Figures
h_x=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
    hold on;
    plot(t*WK.f,des.x_fit{ii}(t), 'k');
    patch_downstroke(h_x,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$x$','interpreter','latex');
% print(h_x, 'hover_control_pos', '-depsc');

h_x_dot=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
    hold on;
    plot(t*WK.f,des.x_dot_fit{ii}(t), 'k');
%     legend('Actual value', 'Fitted ideal value');
    patch_downstroke(h_x_dot,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\dot x$','interpreter','latex');
% print(h_x_dot, 'hover_control_vel', '-depsc');

h_err = figure;
subplot(3,2,1);
plot(t*WK.f, pos_err(1,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
subplot(3,2,3);
plot(t*WK.f, pos_err(2,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
ylabel('$\Delta x$','interpreter','latex');
subplot(3,2,5);
plot(t*WK.f, pos_err(3,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
%
subplot(3,2,2);
plot(t*WK.f, f_a(1,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
subplot(3,2,4);
plot(t*WK.f, f_a(2,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
ylabel('$f_a$','interpreter','latex');
subplot(3,2,6);
plot(t*WK.f, f_a(3,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
xlabel('$t/T$','interpreter','latex');

%%
% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

end

function [X_dot R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau Euler_R Euler_R_dot pos_err]= eom(INSECT, WK_R, WK_L, t, X, des, gains, i, x0)
x=X(1:3);
x_dot=X(4:6);
int_d_x=X(7:9);

% Control design
d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
for j=1:3
%     d_x(j) = des.x_fit{j}(t) - (x(j) - x0(j));
    d_x(j) = des.x_fit{j}(t) - x(j);
    d_x_dot(j) = des.x_dot_fit{j}(t) - x_dot(j);
end
pos_err = INSECT.m*(gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x);
dphi_m = - pos_err(3) / des.df_a_3_by_dphi_m;
dtheta_m = (pos_err(1) + dphi_m * des.df_a_1_by_dphi_m) / des.df_a_1_by_dtheta_m;
dpsi_m = pos_err(2) / des.df_a_2_by_dpsi_m;

WK_R.phi_m = WK_R.phi_m + dphi_m;
WK_L.phi_m = WK_L.phi_m + dphi_m;
WK_R.theta_m = WK_R.theta_m + dtheta_m;
WK_L.theta_m = WK_L.theta_m + dtheta_m;
WK_R.psi_m = WK_R.psi_m + dpsi_m;
WK_L.psi_m = WK_L.psi_m - dpsi_m;

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R Q_L W_R W_L W_R_dot W_L_dot] = wing_attitude(WK_R.beta, Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);

[R W W_dot theta_B] = body_attitude(t, WK_R.f, WK_R); % body
[Q_A W_A W_A_dot theta_A] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen

[L_R L_L D_R D_L M_R M_L ...
    F_rot_R F_rot_L M_rot_R M_rot_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
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
f_a_1=f_a(1:3);
f_a_2=f_a(4:15);

% gravitational force and moment
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:3);
f_g_2=f_g(4:15);

% Euler-Lagrange equation
xi_1=[x_dot];
xi_2=[W; W_R; W_L; W_A];
xi_2_dot=[W_dot; W_R_dot; W_L_dot; W_A_dot];

[JJ KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

[JJ_11 JJ_12 JJ_21 JJ_22] = inertia_sub_decompose_3_12(JJ);
[LL_11 LL_12 LL_21 LL_22] = inertia_sub_decompose_3_12(LL);
[co_ad_11, ~, ~, co_ad_22] = inertia_sub_decompose_3_12(co_ad);

xi_1_dot = JJ_11\( -JJ_12*xi_2_dot -LL_11*xi_1 - LL_12*xi_2 + f_a_1 + f_g_1);

f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1 + JJ_22*xi_2) ...
    + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
f_tau = [zeros(3,1); f_tau_2];
tau = blkdiag(zeros(3), Q_R, Q_L, Q_A)*f_tau_2;

X_dot=[xi_1; xi_1_dot; d_x;];
end
