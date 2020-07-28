function sim_QS_xR_hover_control
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
des = load('sim_QS_xR_hover.mat', 'INSECT', 't', 'x', 'x_dot','R','W',...
    'x0', 'x_dot0', 'R0', 'W0', 'WK');

filename = 'sim_QS_xR_hover_control';
INSECT = des.INSECT;
WK = des.WK;
des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1);
des.R_fit = cell(3, 3); des.W_fit = cell(3, 1);
% 'fourier8', 'cubicinterp'
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', 'fourier8');
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', 'fourier8');
    des.W_fit{i} = fit(des.t, des.W(i, :)', 'fourier8');
    for j=1:3
        des.R_fit{i,j} = fit(des.t, squeeze(des.R(i, j, :)), 'fourier8');
    end
end

% Values obtained from parametric study
load('parametric_study_xR.mat', 'params');
des.params = params;
% Observing parametric study plots
% des.params.mat_aero(1,[3,5,6],:) = 0;
% des.params.mat_aero(2,[1,2,4,6],:) = 0;
% des.params.mat_aero(3,[3,5,6],:) = 0;
% des.params.mat_aero(3,[1,2,4],1) = 1; % A big number?
% des.params.mat_aero(4,[1,2,4,5,6],:) = 0;
% des.params.mat_aero(5,[3,5,6],:) = 0;
% des.params.mat_aero(6,[1,2,4],:) = 0;

% [f_a, M_a] vs [dphi_ms, dtheta_0s, dphi_mk, dphi_0s, dtheta_0k, dpsi_0k]
if rank(des.params.mat_aero(:,:,1), 1e-10) < 6
    error('Parameters are not able to generate aerodynamic forces and moments');
end

N = 1001;
N_period = 3; %15
err_bound = 1e-4; % Convergence criterion is a f(N, N_period, WK)
N_single = round((N-1)/N_period);
T = N_period/WK.f;
t = linspace(0,T,N);
bound_param = 0.1; % Parameter bound; Use 0.25?

des.x_fit_t = zeros(3, N); des.x_dot_fit_t = zeros(3, N);
des.R_fit_t = zeros(3, 3, N); des.W_fit_t = zeros(3, N);
for i=1:3
    des.x_fit_t(i, :) = des.x_fit{i}(t);
    des.x_dot_fit_t(i, :) = des.x_dot_fit{i}(t);
    des.W_fit_t(i, :) = des.W_fit{i}(t);
    for j=1:3
        des.R_fit_t(i, j, :) = des.R_fit{i,j}(t);
    end
end
% des = rmfield(des, {'WK', 'INSECT', 't', 'x', 'x_dot', 'R', 'W', 'f_tau'}); % 'x_fit', 'x_dot_fit'

%% Gains
pol = poly([-7.8 + 19i, -7.8 - 19i, -0.003]);
if ~all(real(roots(pol)) < 0)
    error('The chosen gains are not suitable');
end
% Optimized gs = [427.1529   15.6076  13.4983];
% gains.Kp_pos = pol(3); gains.Kd_pos = pol(2); gains.Ki_pos = pol(4);
gains.Kp_pos = 4; gains.Kd_pos = 2; gains.Ki_pos = 1;
gains.KR = 4; gains.KOm = 2; gains.KI = 1; gains.cI = 1e-1;

%% Single simulation
rng default;
eps = 1e-3;
% dx0 = rand(3,1)*eps;
dx0 = zeros(3, 1);
dx_dot0 = zeros(3, 1);
x0 = des.x0 + dx0;
x_dot0 = des.x_dot0 + dx_dot0;
int_d_x0 = zeros(3, 1);
R0 = des.R0*expmhat(rand(1)*1e-3*[1,0,0]);
W0 = des.W0;% + rand(3,1)*eps;
int_att0 = zeros(3,1);
X0 = [x0; reshape(R0, 9, 1); x_dot0; W0; int_d_x0; int_att0];
wt = 0;

tic;
[err_pos, N_conv, x, x_dot, int_d_x, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, err_xR, dang] =  simulate_control(gains, WK, INSECT, des, X0, N, ...
    N_single, N_period, t, wt, bound_param, err_bound);
toc;
%% Monte Carlo
% N_sims = 10000;
% eps = 3e0;
% wts = [0, 0.1];
% [x_pert, err_pos, N_conv] = monte_carlo(N_sims, eps, wts, gains, WK, ...
%             INSECT, des, N, N_single, N_period, t, bound_param, err_bound);

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

function [x_pert, err_pos, N_conv] = monte_carlo(N_sims, eps, wts, gains, WK, ...
            INSECT, des, N, N_single, N_period, t, bound_param, err_bound)
%%
N_wts = length(wts);
err_pos = zeros(N_sims, N_wts);
N_conv = zeros(N_sims, N_wts);
x_pert = zeros(3, N_sims);
des_x0 = des.x0;
des_x_dot0 = des.x_dot0;
x_dot0 = des_x_dot0;
int_d_x0 = zeros(3, 1);

tic;
par_pool = gcp;
nWorkers = par_pool.NumWorkers;
ticBytes(par_pool);
parfor n=1:nWorkers
    rng(n);
end
pause(1);
parfor i = 1:N_sims
    % % x-z perturbation
    r = rand(1);
    theta = rand(1) * 2*pi;
    dx = [r*cos(theta); 0; r*sin(theta);]*eps;
    % % y perturbation
%     y = rand(1) - 0.5;
%     dx = [0; y; 0;]*eps;
    %
    x_pert(:, i) = dx;
    x0 = des_x0 + dx;
    X0 = [x0; x_dot0; int_d_x0;];
    for l = 1:N_wts
        wt = wts(l);
        [err_pos(i, l), N_conv(i, l)] =  simulate_control(gains, WK, ...
        INSECT, des, X0, N, N_single, N_period, t, wt, bound_param, err_bound);
%         fprintf('%d, %d\n', l, i);
    end
end
tocBytes(par_pool);
toc;

end

function [err_pos, N_conv, x, x_dot, int_d_x, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, err_xR, dang] =  simulate_control(gains, WK, INSECT, des, X0, N, N_single, N_period, t, wt, bound_param, err_bound)
%%
X = zeros(N, 24);
X(1, :) = X0;
dt = t(2) - t(1);
% % Euler Method
% for i=1:(N-1)
%     [X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_A(i) ...
%         W(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
%         W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
%         f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) err_xR(:, i) dang(:, i)]... 
%         = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), wt, bound_param);
% %     theta_B(i) = R2axang(R(:,:,i));
%     X(i+1, 1:3) = X(i, 1:3) + dt * X_dot(1:3, i)';
%     X(i+1, 4:12) = reshape(R(:,:,i)*expmhat(dt*W(:,i)), 9, 1);
%     X(i+1, 13:24) = X(i, 13:24) + dt * X_dot(13:24, i)';
% %     X(i+1, :) = X(i, :) + dt * X_dot(:, i)';
% end
% i = i + 1;
% [X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_A(i) ...
%         W(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
%         W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
%         f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) err_xR(:, i) dang(:, i)]... 
%         = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), wt, bound_param);

% ODE45
[t, X]=ode45(@(t,X) eom_control(INSECT, WK, WK, t, X, des, gains, 0, X0(1:3), wt, bound_param), ...
    t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
for k=1:N    
    [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),...
        theta_A(k), W(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
        W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
        M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k), ...
        Euler_R(:,k), Euler_R_dot(:,k), err_xR(:, k), dang(:, k)]= ...
        eom_control(INSECT, WK, WK, t(k), X(k,:)', des, gains, k, X0(1:3), wt, bound_param);
end

t_single = (max(t)/N_period);
err = vecnorm(err_xR, 2, 1);
N_conv = Inf;
for j = 1:N_single:(N-N_single)
    err_pos = sum(err(j:(j+N_single))) * dt / t_single;
    if err_pos <= err_bound
        N_conv = (j-1)/N_single + 1;
        break;
    end
end
err_pos = sum(err(:, N-N_single:N)) * dt / t_single;
% err_pos = norm(pos_err(:, end));

x=X(:,1:3)';
x_dot=X(:,13:15)';
int_d_x=X(:,19:21)';

end

function [X_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, err_xR, dang]= eom_control(INSECT, WK_R, WK_L, t, X, des, gains, i, x0, wt, bound_param)
%% Dynamics along with the designed control
x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);
int_d_x=X(19:21);
int_att=X(22:24);

%% Ideal values
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
[JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_6_9(JJ);
C=[zeros(3,9);
    -Q_R -Q_L -Q_A];
% [JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(4:6), W, W_A);
% f_abd = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
% % f_abd(isnan(f_abd)) = 0;
% 
[L_R, L_L, D_R, D_L, M_R, M_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
F_R = L_R + D_R;
F_L = L_L + D_L;
f_a=[R*Q_R*F_R + R*Q_L*F_L;
    hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L + Q_R*M_R + Q_L*M_L];
% % f_a(isnan(f_a)) = 0;
% f_total = f_a + f_abd; % Use this for gettings signs

%% Control design
d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
des_R = zeros(3, 3); des_W = zeros(3, 1);
for k=1:3
    d_x(k) = des.x_fit{k}(t) - x(k);
    d_x_dot(k) = des.x_dot_fit{k}(t) - x_dot(k);
    des_W(k) = des.W_fit{k}(t);
    for j=1:3
        des_R(k, j) = des.R_fit{k,j}(t);
    end
end
e_R = 0.5*vee(des_R'*R - R'*des_R);
e_Om = W - R'*des_R*des_W;

% d_x = des.x(:, i) - x;
% d_x_dot = des.x_dot(:, i) - x_dot;
% e_R = 0.5*vee(des.R(:, :, i)'*R - R'*des.R(:, :, i));
% e_Om = W - R'*des.R(:, :, i)*des.W(:,i);

% d_x = des.x_fit_t(:, i) - x;
% d_x_dot = des.x_dot_fit_t(:, i) - x_dot;
% e_R = 0.5*vee(des.R_fit_t(:, :, i)'*R - R'*des.R_fit_t(:, :, i));
% e_Om = W - R'*des.R_fit_t(:, :, i)*des.W_fit_t(:,i);
int_att_dot = e_Om + gains.cI * e_R;

pos_err = (gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x); % INSECT.m*()
att_err = (- gains.KR*e_R - gains.KOm*e_Om - gains.KI*int_att);
JJ_xR = (JJ_11 - C*JJ_21);
err_xR = JJ_xR*[pos_err; att_err];
% err_xR = [JJ_xR(1:3,1:3)*pos_err; JJ_xR(4:6,4:6)*att_err];
% err_xR = JJ_11 * [pos_err; att_err];
% err_xR = [JJ_11(1:3,1:3)*pos_err; JJ_11(4:6,4:6)*att_err];

sign_f_a = round((3 - sign(f_a))/2);
params = des.params.mat_aero(1:6, 1:6, 1);
for k=1:6
    params(k,:) = des.params.mat_aero(k,:,sign_f_a(k));
end

dang = zeros(7,1);
if wt == 0
    dang(1:6) = params \ err_xR;
else
    % Minimum norm solution
    temp_A = temp_A(1:2,:);
%     pos_err = tanh(pos_err) * 1e-2; %(norm(temp_A, 'fro') * bound_param);
%     rhs = [pos_err(1); pos_err(3); 0];
    dang = temp_A' * ((temp_A*temp_A') \ rhs(1:2));
end

idx = abs(dang) > bound_param;
dang(idx) = bound_param * sign(dang(idx));

% dang(:) = 0;
% dang(1) = 0.1;%*sin(t*WK_R.f*2*pi);

WK_R.phi_m = WK_R.phi_m + dang(1) + dang(3);
WK_L.phi_m = WK_L.phi_m + dang(1) - dang(3);
WK_R.phi_0 = WK_R.phi_0 + dang(4);
WK_L.phi_0 = WK_L.phi_0 + dang(4);
WK_R.theta_0 = WK_R.theta_0 + dang(2) + dang(5);
WK_L.theta_0 = WK_L.theta_0 + dang(2) - dang(5);
WK_R.psi_m = WK_R.psi_m + dang(6);
WK_L.psi_m = WK_L.psi_m - dang(6);
WK_R.theta_A_m = WK_R.theta_A_m + dang(7);
WK_L.theta_A_m = WK_L.theta_A_m + dang(7);

X = X(1:18);
[X_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau, Euler_R, Euler_R_dot] = eom_QS_xR(INSECT, WK_R, WK_L, t, X);

% X = X(1:18);
% u_control = err_xR;
% [X_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, ...
%     W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
%     f_tau, tau, Euler_R, Euler_R_dot] = eom_QS_xR_ideal(INSECT, WK_R, WK_L, t, X, u_control);

X_dot=[X_dot(1:18); d_x; int_att_dot];

end
