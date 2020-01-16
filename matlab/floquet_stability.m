function floquet_stability
% Study the stability of a periodic trajecoty using floquet theory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data');

%% Linearized dynamics
filename=append('sim_QS_x_hover_','stability_data');
load('sim_QS_x_hover.mat');
INSECT.scale=1e-2;
INSECT.name='MONARCH';
WK.ab_type='varying';
WK.bo_type='varying';
% load('sim_QS_x_hover_hawkmoth.mat');
% INSECT.name='NOT_MONARCH';

N_sims = 10;
conv_rate_osc = zeros(N_sims, 6);
var_name_to_save = 'conv_rate_osc';

N=1001;
N_periods=10;
T=N_periods/WK.f;
ix_d = (N-1)/N_periods;
t=linspace(0,T,N);
dt = T/(N-1);
epsilon = 1e0;

[delta_mat, F_linear] = sim_perturbation(INSECT, WK, X0, N, t, epsilon);

B = zeros(6, 6, 1+ix_d);
start_ix = max(1, round((N_periods-2)/N_periods * N));
for i=start_ix:(start_ix+ix_d)
    j = i-start_ix+1;
%         d_F = (F_linear(:, :, k) - F_linear(:, :, k+ix_d))./ F_linear(:, :, k);
%         d_F(max(abs(F_linear(:, :, k+ix_d)), abs(F_linear(:, :, k))) < 1e-10) = 0;
%         if(~all(abs(d_F) < 1e-2, [1, 2]))
%             disp(d_F)
%         end
    B(:, :, j) = delta_mat(:, :, i) \ delta_mat(:, :, i+ix_d);
end
B = B(:, :, round(end/2));
[e_vecs, rhos] = eig(B);
mus = log(diag(rhos)) * WK.f;

char_soln_mat = zeros(6, 6, N);
per_val_mat = zeros(6, 6, N);

for i=1:N
    char_soln_mat(:, :, i) = delta_mat(:, :, i) * e_vecs;
    Y_0 = diag(exp(mus*t(i)));
    per_val_mat(:, :, i) = char_soln_mat(:, :, i) * inv(Y_0);
end

omega = 2*pi*WK.f;
stop_idx = N;

for c_ix=4:6 % Column index for perturbation direction
%     delta_g_mag = vecnorm(reshape(delta_mat(1:3, c_ix, :), 3, N));
%     delta_xi_mag = vecnorm(reshape(delta_mat(4:6, c_ix, :), 3, N));
    conv_rate_osc(:, c_ix) = repmat(mus(c_ix), [N_sims, 1]);
end

% save('sim_QS_x_hover_conv_rate', var_name_to_save, '-append');

%% Figures
% time=t*WK.f;
% figure;
% xlabel('$t/T$','interpreter','latex');
% subplot(2, 1, 1);
% plot(time(1:stop_idx), delta_g_mag(1:stop_idx));
% ylabel('$\delta x$','interpreter','latex');
% hold on;
% subplot(2, 1, 2);
% plot(time(1:stop_idx), delta_xi_mag(1:stop_idx));
% ylabel('$\delta \dot{x}$','interpreter','latex');
% % print('sim_QS_x_hover_stability', '-depsc');

% h_floq = figure;
% h_floq.PaperUnits = 'inches';
% h_floq.PaperPosition = [0 0 11 11];
% for i=1:3
%     subplot(3, 2, 2*i-1);
%     plot(t*WK.f, squeeze(char_soln_mat(:, i+3, :)));
%     lgd = legend(string(1:6));
% %     title(lgd, 'Component index');
%     l = "$\mathbf{x}_" + string(i+3) + "(t)$" + " with " + "$\mu_" + string(i+3) + " = " + string(mus(i+3)) + "$";
%     ylabel(l, 'interpreter', 'latex');
%     xlabel('$t/T$', 'interpreter', 'latex');
% 
%     j = 2*i;
%     subplot(3, 2, j);
%     plot(t*WK.f, squeeze(per_val_mat(:, i+3, :)));
%     lgd = legend(string(1:6));
%     l = "$\mathbf{p}_" + string(i+3) + "(t)$" + " with " + "$\mu_" + string(i+3) + " = " + string(mus(i+3)) + "$";
%     ylabel(l, 'interpreter', 'latex');
%     xlabel('$t/T$', 'interpreter', 'latex');
% end
% print(h_floq, 'hover_char_soln', '-depsc', '-r0');

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

%%
function [delta_mat, F_linear] = sim_perturbation(INSECT, WK, X0, N, t, epsilon)

delta0 = diag(rand(6, 1))/epsilon;
% delta0 = rand(6, 6)/epsilon;
X0 = [X0(1:6); reshape(delta0, 36, 1);];

[t X]=ode45(@(t,X) eom(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

F_linear=zeros(6, 6, N);
for k=1:N
    [~, F_linear(:, :, k)] = eom(INSECT, WK, WK, t(k), X(k, :)');
end

x=X(:,1:3)';
x_dot=X(:,4:6)';

% delta=X(:,7:12)';
% delta_g_mag = sqrt(diag(delta(1:3, :)'*delta(1:3, :)));
% delta_xi_mag = sqrt(diag(delta(4:6, :)'*delta(4:6, :)));
delta_mat=reshape(X(:,7:42)', 6, 6, N);

end

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom(INSECT, WK_R, WK_L, t, X)
x=X(1:3);
x_dot=X(4:6);
% delta=X(7:12);
delta_mat=reshape(X(7:42), 6, 6);

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R Q_L W_R W_L W_R_dot W_L_dot] = wing_attitude(WK_R.beta, Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);

% [R W W_dot theta_B] = body_attitude(t,WK_R.f); %time-varying thorax
% [Q_A W_A W_A_dot theta_A] = abdomen_attitude(t,WK_R.f); % abdomen

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

F_linear = zeros(6, 6);
[d_L_R d_L_L d_D_R d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
d_F_R=d_L_R+d_D_R;
d_F_L=d_L_L+d_D_L;
F_linear(1:3, 4:6) = eye(3);
F_linear(4:6, 4:6) = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L);
X_dot=[xi_1; xi_1_dot; reshape(F_linear*delta_mat, 36, 1);];

end

