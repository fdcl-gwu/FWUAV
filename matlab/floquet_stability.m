evalin('base','clear all');

%% Linearized dynamics
load('sim_QS_x_hover_multistart_vary_vel_and_freq_final_better_no_tol_correct_abdomen_osc_less_freq.mat');
INSECT.scale=1e-2;
INSECT.name='MONARCH';
% load('sim_QS_x_hover_hawkmoth.mat');
% INSECT.name='NOT_MONARCH';

N_sims = 10;
conv_rate_osc = zeros(N_sims, 6);
var_name_to_save = 'conv_rate_osc';

N=1001; % 3001
N_periods=10;
T=N_periods/WK.f;
% find(t==T/2);
ix_d = (N-1)/N_periods;
t=linspace(0,T,N);
dt = T/(N-1);
epsilon = 1e0;

% poolobj = parpool(4);
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

% polarplot(complex(eig_vals_B(:, end)), 'r*');
% print('char_multipliers', '-depsc');

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
%     ft=fittype('a*exp(-b*x)');
%     fit_value=fit(t(1:stop_idx)', delta_xi_mag(1:stop_idx)', ft, 'StartPoint', [1, 0.1]);
%     fit_delta_xi_mag = fit_value.a * exp(-fit_value.b * t);
%     conv_rate = fit_value.b;
%     conv_rate_osc(:, c_ix) = repmat(conv_rate, [N_sims, 1]);
    conv_rate_osc(:, c_ix) = repmat(mus(c_ix), [N_sims, 1]);
end

% delete(poolobj)

% save('sim_QS_x_hover_conv_rate', var_name_to_save, '-append');

time=t*WK.f;

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

% filename='sim_QS_x_hover';
% filename=append('sim_QS_x_hover_','stability_data');

% % Get a list of all variables
% allvars = whos;
% % Identify the variables that ARE NOT graphics handles. This uses a regular
% % expression on the class of each variable to check if it's a graphics object
% tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% % Pass these variable names to save
% save(filename, allvars(tosave).name)
% evalin('base',['load ' filename]);

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

function fourier_analysis()
% c_ix = 6;
% per_val = vecnorm(squeeze(per_val_mat(:, c_ix, :)));

% omega = 2*pi*WK.f;
% ft_per=fittype('a*sin(omega*x + b) + c');
% fit_val_per=fit(t(1:stop_idx)', per_val(1:stop_idx)', ft_per, ...
%     'StartPoint', [1, 1, 1, 1], 'Lower', [-Inf, -Inf, -Inf, omega], 'Upper', [Inf, Inf, Inf, omega]);
% fit_per_val = fit_val_per.a * sin(omega*t(1:stop_idx) + fit_val_per.b) + fit_val_per.c;
% 
% per_fft = fft(per_val-fit_val_per.c);
% P2 = abs(per_fft/N);
% P1 = P2(1:round(N/2));
% P1(2:end-1) = 2*P1(2:end-1);
% f = WK.f*(1:round(N/2))/N;
% figure
% plot(f, P1);
% axis auto;
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

[R W W_dot theta_B] = body_attitude(WK_R.theta_B); % body
[Q_A W_A W_A_dot theta_A] = abdomen_attitude(t, WK_R, 'designed'); % abdomen

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

function [JJ_11 JJ_12 JJ_21 JJ_22] = inertia_sub_decompose_3_12(JJ)
JJ_11 = JJ(1:3,1:3);
JJ_12 = JJ(1:3,4:15);
JJ_21 = JJ(4:15,1:3);
JJ_22 = JJ(4:15,4:15);
end
    
function [JJ KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
[JJ_R KK_R] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R, Q_R, x_dot, W, W_R);
[JJ_L KK_L] = inertia_wing_sub(INSECT.m_L, INSECT.mu_L, INSECT.nu_L, INSECT.J_L, R, Q_L, x_dot, W, W_L);
[JJ_A KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, x_dot, W, W_A);

JJ=zeros(15,15);
JJ(1:3,1:3) = INSECT.m_B*eye(3) + JJ_R(1:3,1:3) + JJ_L(1:3,1:3) + + JJ_A(1:3,1:3);
JJ(1:3,4:6) = JJ_R(1:3,4:6) + JJ_L(1:3,4:6) + JJ_A(1:3,4:6);
JJ(1:3,7:9) = JJ_R(1:3,7:9);
JJ(1:3,10:12) = JJ_L(1:3,7:9);
JJ(1:3,13:15) = JJ_A(1:3,7:9);

JJ(4:6,1:3) = JJ(1:3,4:6)';
JJ(4:6,4:6) = INSECT.J_B + JJ_R(4:6,4:6) + JJ_L(4:6,4:6) + + JJ_A(4:6,4:6);
JJ(4:6,7:9) = JJ_R(4:6,7:9);
JJ(4:6,10:12) = JJ_L(4:6,7:9);
JJ(4:6,13:15) = JJ_A(4:6,7:9);

JJ(7:9,1:3) = JJ(1:3,7:9)';
JJ(7:9,4:6) = JJ(4:6,7:9)';
JJ(7:9,7:9) = JJ_R(7:9,7:9);

JJ(10:12,1:3) = JJ(1:3,10:12)';
JJ(10:12,4:6) = JJ(4:6,10:12)';
JJ(10:12,10:12) = JJ_L(7:9,7:9);

JJ(13:15,1:3) = JJ(1:3,13:15)';
JJ(13:15,4:6) = JJ(4:6,13:15)';
JJ(13:15,13:15) = JJ_A(7:9,7:9);

KK=zeros(15,15);
KK(1:3,4:6) = KK_R(1:3,4:6) + KK_L(1:3,4:6) + KK_A(1:3,4:6);
KK(1:3,7:9) = KK_R(1:3,7:9);
KK(1:3,10:12) = KK_L(1:3,7:9);
KK(1:3,13:15) = KK_A(1:3,7:9);

KK(4:6,4:6) = KK_R(4:6,4:6) + KK_L(4:6,4:6) + KK_A(4:6,4:6);
KK(4:6,7:9) = KK_R(4:6,7:9);
KK(4:6,10:12) = KK_L(4:6,7:9);
KK(4:6,13:15) = KK_A(4:6,7:9);

KK(7:9,4:6) = KK_R(7:9,4:6);
KK(7:9,7:9) = KK_R(7:9,7:9);

KK(10:12,4:6) = KK_L(7:9,4:6);
KK(10:12,10:12) = KK_L(7:9,7:9);

KK(13:15,4:6) = KK_A(7:9,4:6);
KK(13:15,13:15) = KK_A(7:9,7:9);
end

function [JJ KK] = inertia_wing_sub(m, mu, nu, J, R, Q, x_dot, W, W_i)
R_dot=R*hat(W);
Q_dot=Q*hat(W_i);

JJ=zeros(9,9);

JJ(1:3,1:3)=m*eye(3);
JJ(1:3,4:6)=-m*R*(hat(mu)+hat(Q*nu));
JJ(1:3,7:9)=-m*R*Q*hat(nu);

JJ(4:6,1:3)=JJ(1:3,4:6)';
JJ(4:6,4:6)=m*hat(mu)'*hat(mu)+Q*J*Q'+m*(hat(mu)'*hat(Q*nu)+hat(Q*nu)'*hat(mu));
JJ(4:6,7:9)=Q*J+m*hat(mu)'*Q*hat(nu);

JJ(7:9,1:3)=JJ(1:3,7:9)';
JJ(7:9,4:6)=JJ(4:6,7:9)';
JJ(7:9,7:9)=J;

KK=zeros(9,9);

KK(1:3,4:6) = m*R*hat((hat(mu)+hat(Q*nu))*W) + m*R*hat(Q*hat(nu)*W_i);
KK(1:3,7:9) = -m*R*hat(W)*Q*hat(nu) + m*R*Q*hat(hat(nu)*W_i);
KK(4:6,4:6) = m*(hat(mu)+hat(Q*nu))*hat(R'*x_dot);
KK(4:6,7:9) = m*hat(R'*x_dot)*Q*hat(nu) - Q*hat(J*Q'*W) + Q*J*hat(Q'*W) ...
    -m*hat(mu)*hat(W)*Q*hat(nu) - m* hat(hat(mu)*W)*Q*hat(nu) ...
    -Q*hat(J*W_i) + m*hat(mu)*Q*hat(hat(nu)*W_i);
KK(7:9,4:6) = m*hat(nu)*Q'*hat(R'*x_dot);
KK(7:9,7:9) = m*hat(nu)*hat(Q'*R'*x_dot) + J*hat(Q'*W) - m*hat(nu)*hat(Q'*hat(mu)*W);
end

function [U dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A)
e3=[0 0 1]';

mg_B=INSECT.m_B*INSECT.g;
mg_R=INSECT.m_R*INSECT.g;
mg_L=INSECT.m_L*INSECT.g;
mg_A=INSECT.m_A*INSECT.g;

tmp_R = INSECT.mu_R + Q_R*INSECT.nu_R;
tmp_L = INSECT.mu_L + Q_L*INSECT.nu_L;
tmp_A = INSECT.mu_A + Q_A*INSECT.nu_A;

U_B = -mg_B*e3'*x;
U_R = -mg_R*e3' * (x + R*tmp_R);
U_L = -mg_L*e3' * (x + R*tmp_L);
U_A = -mg_A*e3' * (x + R*tmp_A);
U = U_B + U_R + U_L + U_A;

dU = [-(INSECT.m_B + INSECT.m_R + INSECT.m_L + INSECT.m_A) * INSECT.g * e3;
    mg_R*hat(R'*e3)*tmp_R + mg_L*hat(R'*e3)*tmp_L + mg_A*hat(R'*e3)*tmp_A;
    mg_R*hat(Q_R'*R'*e3)*INSECT.nu_R;
    mg_L*hat(Q_L'*R'*e3)*INSECT.nu_L;
    mg_A*hat(Q_A'*R'*e3)*INSECT.nu_A];
end
