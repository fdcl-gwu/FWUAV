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
N_study = 7;
eps = linspace(-0.1, 0.1, N_params);
f_aero = zeros(N_study, 3, N, N_params);
M_aero = zeros(N_study, 3, N, N_params);
f_a_m = zeros(N_study, 3, N_params, 2);
M_a_m = zeros(N_study, 3, N_params, 2);

[t,X]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
for k=1:N
    [~, R(:, :, k), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, F_R(:, k), F_L(:, k), ~, ~, f_a(:,k)]= eom_QS_xR(INSECT, WK, WK, t(k), X(k,:)');
end

parfor i=1:N_params
    [f_aero(:,:,:,i), M_aero(:,:,:,i)] = param_study(INSECT, WK, eps, i, t, X0, N, R, F_R);
    f_aero1 = f_aero(:, :, :, i); f_aero2 = f_aero(:, :, :, i);
    f_aero1(f_aero1 < 0) = 0;
    f_aero2(f_aero2 > 0) = 0;
    
    f_a_m_temp = zeros(size(f_a_m(:, :, i, :)));
    f_a_m_temp(:, :, 1, 1) = mean(f_aero1, 3);
%     f_a_m_temp(:, :, 1, 1) = mean(f_aero(:,:,:,i), 3);
    f_a_m_temp(:, :, 1, 2) = mean(f_aero2, 3);
    f_a_m(:, :, i, :) = f_a_m_temp;

%     M_a_m(:, :, i, 1) = mean(M_aero(:,:,:,i), 3);
    M_aero1 = M_aero(:, :, :, i); M_aero2 = M_aero(:, :, :, i);
    M_aero1(M_aero1 < 0) = 0;
    M_aero2(M_aero2 > 0) = 0;
    M_a_m_temp = zeros(size(M_a_m(:, :, i, :)));
    M_a_m_temp(:, :, 1, 1) = mean(M_aero1, 3);
    M_a_m_temp(:, :, 1, 2) = mean(M_aero2, 3);
    M_a_m(:, :, i, :) = M_a_m_temp;
end
% f_a_m = abs(f_a_m);

params.df_a_1_by_dphi_ms = [lin_reg(eps', squeeze(f_a_m(1,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(1,1,:,2)))];
params.df_a_3_by_dphi_ms = lin_reg(eps', squeeze(f_a_m(1,3,:,2)));
params.df_a_1_by_dtheta_0s = [lin_reg(eps', squeeze(f_a_m(2,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(2,1,:,2)))];
params.df_a_3_by_dtheta_0s = lin_reg(eps', squeeze(f_a_m(2,3,:,2)));
params.df_a_2_by_dphi_mk = max_abs([lin_reg(eps(eps'>0)', squeeze(f_a_m(3,2,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(f_a_m(3,2,eps'<0,1)))]);
% But which one to use is undecided - use the one with the larger magnitude;
params.df_a_1_by_dphi_0s = [lin_reg(eps', squeeze(f_a_m(4,1,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(4,1,:,2)))];
params.df_a_3_by_dphi_0s = lin_reg(eps', squeeze(f_a_m(4,3,:,2)));
params.df_a_2_by_dtheta_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(f_a_m(5,2,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(f_a_m(5,2,eps'<0,1)))]);
params.df_a_2_by_dpsi_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(f_a_m(6,2,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(f_a_m(6,2,eps'<0,1)))]);
%
% params.df_a_1_by_dtheta_A_m = [lin_reg(eps', squeeze(f_a_m(4,1,:,1))), ...
%     lin_reg(eps', squeeze(f_a_m(4,1,:,2)))];
% params.df_a_3_by_dtheta_A_m = [lin_reg(eps', squeeze(f_a_m(4,3,:,1))), ...
%     lin_reg(eps', squeeze(f_a_m(4,3,:,2)))];
%
params.dM_a_2_by_dphi_ms = [lin_reg(eps', squeeze(M_a_m(1,2,:,1))), ...
    lin_reg(eps', squeeze(f_a_m(1,2,:,2)))];
params.dM_a_2_by_dtheta_0s = [lin_reg(eps', squeeze(M_a_m(2,2,:,1))), ...
    lin_reg(eps', squeeze(M_a_m(2,2,:,2)))];
params.dM_a_1_by_dphi_mk = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(3,1,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(3,1,eps'<0,1)))]);
params.dM_a_3_by_dphi_mk = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(3,3,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(3,3,eps'<0,1)))]);
params.dM_a_2_by_dphi_0s = [lin_reg(eps', squeeze(M_a_m(4,2,:,1))), ...
    lin_reg(eps', squeeze(M_a_m(4,2,:,2)))];
params.dM_a_1_by_dtheta_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(5,1,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(5,1,eps'<0,1)))]);
params.dM_a_3_by_dtheta_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(5,3,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(5,3,eps'<0,1)))]);
params.dM_a_1_by_dpsi_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(6,1,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(6,1,eps'<0,1)))]);
params.dM_a_3_by_dpsi_0k = max_abs([lin_reg(eps(eps'>0)', squeeze(M_a_m(6,3,eps'>0,1))), ...
    lin_reg(eps(eps'<0)', squeeze(M_a_m(6,3,eps'<0,1)))]);
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

function [f_aero, M_aero] = param_study(INSECT, WK, eps, i, t, X0, N, R, F_R)
%%
    f_aero = zeros(7, 3, N);
    M_aero = zeros(7, 3, N);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.phi_m = WK_R.phi_m + eps(i);
    WK_L.phi_m = WK_L.phi_m + eps(i);
    [f_aero(1, :, :), M_aero(1, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.theta_0 = WK_R.theta_0 + eps(i);
    WK_L.theta_0 = WK_L.theta_0 + eps(i);
    [f_aero(2, :, :), M_aero(2, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.phi_m = WK_R.phi_m + eps(i);
    WK_L.phi_m = WK_L.phi_m - eps(i);
    [f_aero(3, :, :), M_aero(3, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.phi_0 = WK_R.phi_0 + eps(i);
    WK_L.phi_0 = WK_L.phi_0 + eps(i);
    [f_aero(4, :, :), M_aero(4, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.theta_0 = WK_R.theta_0 + eps(i);
    WK_L.theta_0 = WK_L.theta_0 - eps(i);
    [f_aero(5, :, :), M_aero(5, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
%     WK_R.phi_0 = WK_R.phi_0 + eps(i);
%     WK_L.phi_0 = WK_L.phi_0 - eps(i);
    WK_R.psi_0 = WK_R.psi_0 + eps(i);
    WK_L.psi_0 = WK_L.psi_0 - eps(i);
    [f_aero(6, :, :), M_aero(6, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
    %
    WK_R = WK;  WK_L = WK;
    WK_R.theta_A_m = WK_R.theta_A_m + eps(i);
    WK_L.theta_A_m = WK_L.theta_A_m + eps(i);
    [~, M_aero(7, :, :), f_aero(7, :, :)] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R, F_R);
end

function [f_aero, M_aero, f_abd] = aerodynamic_param(INSECT, WK_R, WK_L, t, X0, N, R_id, F_R_id)
%%
    [t,X]=ode45(@(t,X) eom_QS_xR(INSECT, WK_R, WK_L, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    f_aero = zeros(3, N); f_abd = zeros(3, N); M_aero = zeros(3, N);
    for k=1:N    
        [~, R, Q_R, Q_L, Q_A, ~, W, ~, ~, ~, ~, W_A, W_A_dot, ~, ~, M_R, M_L, f_a, ~, f_tau]= eom_QS_xR(INSECT, WK_R, WK_L, t(k), X(k,:)');
        [JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(k, 4:6)', W, W_A);
        f_abd(:, k) = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
        f_aero(:, k) = f_a(1:3);
%         M_aero(:, k) = f_a(4:6);
%         M_aero(:, k) = f_a(4:6) + f_tau(4:6);
        M_aero(:, k) = f_a(4:6) + Q_R*M_R + Q_L*M_L;
    end
    f_total = f_aero + f_abd;
end

function out = max_abs(inp)
[~, I] = max(abs(inp));
out = inp(I);
end
