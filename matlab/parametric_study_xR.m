function parametric_study_xR
% Simulate the effects on aerodynamic forces when changing wing kinematic
% parameters considering (x,R) dynamics

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
load('sim_QS_xR_hover.mat', 'INSECT', 'WK', 'X0');
filename='parametric_study_xR';

N=501;
T=1/WK.f;
t=linspace(0,T,N);

N_params = 50;
N_study = 4;
eps = linspace(-0.1, 0.1, N_params);
f_aero = zeros(N_study, 3, N, N_params);
f_a_m = zeros(N_study, 3, N_params, 2);

[t,X]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
for k=1:N
    [~, R(:, :, k), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, F_R(:, k), F_L(:, k), ~, ~, f_a(:,k)]= eom_QS_xR(INSECT, WK, WK, t(k), X(k,:)');
end

parfor i=1:N_params
    f_aero(:, :, :, i) = param_study(INSECT, WK, eps, i, t, X0, N, R, F_R);
    f_aero1 = f_aero(:, :, :, i); f_aero2 = f_aero(:, :, :, i);
    f_aero1(f_aero1 < 0) = 0;
    f_aero2(f_aero2 > 0) = 0;
    
    f_a_m_temp = zeros(size(f_a_m(:, :, i, :)));
    f_a_m_temp(:, :, 1, 1) = mean(f_aero1, 3);
    f_a_m_temp(:, :, 1, 2) = mean(f_aero2, 3);
    f_a_m(:, :, i, :) = f_a_m_temp;
end
% f_a_m = abs(f_a_m);

params.df_a_1_by_dphi_m = [lin_reg(eps', squeeze(f_a_m(1,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(1,1,:,2)))]; % dphi_m_R > 0, dphi_m_L > 0
params.df_a_3_by_dphi_m = lin_reg(eps', squeeze(f_a_m(1,3,:,2))); % dphi_m_R > 0, dphi_m_L > 0
params.df_a_1_by_dtheta_0 = [lin_reg(eps', squeeze(f_a_m(2,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(2,1,:,2)))]; % dtheta_0 > 0
params.df_a_3_by_dtheta_0 = lin_reg(eps', squeeze(f_a_m(2,3,:,2))); % dtheta_0 > 0
params.df_a_2_by_dphi_m = [lin_reg(eps(eps'>0)', squeeze(f_a_m(3,2,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(f_a_m(3,2,eps'<0,1)))]; % dphi_m_R > 0, dphi_m_L < 0;
% But which one to use is undecided - depends on the stroke probably;
% downstroke : use +ve value, upstroke : use -ve value
params.df_a_1_by_dtheta_A_m = [lin_reg(eps', squeeze(f_a_m(4,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(4,1,:,2)))]; % dtheta_A_m > 0
params.df_a_3_by_dtheta_A_m = [lin_reg(eps', squeeze(f_a_m(4,3,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(4,3,:,2)))]; % dtheta_A_m > 0

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

function f_aero = param_study(INSECT, WK, eps, i, t, X0, N, R, F_R)
%%
    f_aero = zeros(4, 3, N);
%     %
%     WK_R = WK;  WK_L = WK;
%     WK_R.phi_m = WK_R.phi_m + eps(i);
%     WK_L.phi_m = WK_L.phi_m + eps(i);
%     f_aero(1, :, :) = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
%     %
%     WK_R = WK;  WK_L = WK;
%     WK_R.theta_0 = WK_R.theta_0 + eps(i);
%     WK_L.theta_0 = WK_L.theta_0 + eps(i);
%     f_aero(2, :, :) = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
%     %
%     WK_R = WK;  WK_L = WK;
% %     WK_R.psi_m = WK_R.psi_m + eps(i);
% %     WK_L.psi_m = WK_L.psi_m - eps(i);
%     WK_R.phi_m = WK_R.phi_m + eps(i);
%     WK_L.phi_m = WK_L.phi_m - eps(i);
%     f_aero(3, :, :) = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
%     %
%     WK_R = WK;  WK_L = WK;
%     WK_R.theta_A_m = WK_R.theta_A_m + eps(i);
%     WK_L.theta_A_m = WK_L.theta_A_m + eps(i);
%     [~, f_aero(4, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.phi_m = WK_R.phi_m + eps(i);
    WK_L.phi_m = WK_L.phi_m - eps(i);
    [~, ~, f_aero(1, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.theta_0 = WK_R.theta_0 + eps(i);
    WK_L.theta_0 = WK_L.theta_0 - eps(i);
    [~, ~, f_aero(2, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.psi_m = WK_R.psi_m + eps(i);
    WK_L.psi_m = WK_L.psi_m - eps(i);
    [~, ~, f_aero(3, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.theta_A_m = WK_R.theta_A_m + eps(i);
    WK_L.theta_A_m = WK_L.theta_A_m + eps(i);
    [~, ~, f_aero(4, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
end

function [f_aero, f_abd, M_aero] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R_id, F_R_id)
%%
    [t,X]=ode45(@(t,X) eom_QS_xR(INSECT, WK_R, WK_L, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    f_a = zeros(15, N); f_abd = zeros(3, N); M_aero = zeros(3, N);
    for k=1:N    
        [~, R, Q_R, Q_L, Q_A, ~, W, ~, ~, ~, ~, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a(:,k)]= eom_QS_xR(INSECT, WK_R, WK_L, t(k), X(k,:)');
        [JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(k, 4:6)', W, W_A);
        f_abd(:, k) = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
        M_aero(:, k) = hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    end
    f_aero = f_a(1:3, :);
    f_total = f_aero + f_abd;
end
