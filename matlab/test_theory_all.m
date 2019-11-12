evalin('base','clear all');
load('morp_MONARCH');
INSECT=MONARCH;
load('sim_QS_x_hover_surrogate_zerovel_inter_low_pattern_fmin_accurate_multistart.mat')

%% Linearized dynamics

N=1001; % 3001
T=2/WK.f;
% find(t==T/2);
ix_d = (N-1)/2;
t=linspace(0,T,N);
dt = T/(N-1);
epsilon = 1e2;
% delta0 = zeros(30, 1);
% delta0(4) = rand(1, 1)/epsilon;
delta0 = diag(rand(1, 30))/epsilon;
X0 = [X0; reshape(delta0, 900, 1);];
[t X]=ode45(@(t,X) eom(INSECT, WK, WK, t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,4:6)';
% delta=X(:,7:36)';
delta_mat=reshape(X(:,7:906)', 30, 30, N);
B = zeros(30, 30, 1+ix_d);
for i=1:(1+ix_d)
    B(:, :, i) = delta_mat(:, :, i) \ delta_mat(:, :, i+ix_d);
%     if(rank(delta_mat(:, :, i)) < 30)
%         disp(delta_mat(:, :, i))
%     end
end
time=t*WK.f;
% delta_g_mag = sqrt(diag(delta(1:15, :)'*delta(1:15, :)));
% delta_xi_mag = sqrt(diag(delta(16:30, :)'*delta(16:30, :)));
c_ix = 4;
delta_g_mag = sqrt(diag(reshape(delta_mat(1:15, c_ix, :), 15, N)'*reshape(delta_mat(1:15, c_ix, :), 15, N)));
delta_xi_mag = sqrt(diag(reshape(delta_mat(16:30, c_ix, :), 15, N)'*reshape(delta_mat(16:30, c_ix, :), 15, N)));
stop_idx = N;

figure;
xlabel('$t/T$','interpreter','latex');
subplot(2, 1, 1);
plot(time(1:stop_idx), delta_g_mag(1:stop_idx));
ylabel('$\delta g$','interpreter','latex');
hold on;
subplot(2, 1, 2);
plot(time(1:stop_idx), delta_xi_mag(1:stop_idx));
ylabel('$\delta \xi$','interpreter','latex');

function [X_dot R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom(INSECT, WK_R, WK_L, t, X)
x=X(1:3);
x_dot=X(4:6);
% delta=X(7:36);
delta_mat=reshape(X(7:906), 30, 30);

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R Q_L W_R W_L W_R_dot W_L_dot] = wing_attitude(WK_R.beta, Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);

% [R W W_dot theta_B] = body_attitude(t,WK_R.f); %time-varying thorax
% [Q_A W_A W_A_dot theta_A] = abdomen_attitude(t,WK_R.f); %time-varying abdomen

[R W W_dot theta_B] = body_attitude(15.65*pi/180); % fixed body
[Q_A W_A W_A_dot theta_A] = abdomen_attitude(17.32*pi/180); % fixed abdomen

% [R W W_dot theta_B] = body_attitude(t, WK_R, 'designed'); % body
% [Q_A W_A W_A_dot theta_A] = abdomen_attitude(WK_R.theta_A); %time-varying abdomen

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

F_linear = zeros(30, 30);
% trace_integ = 0;
x_ddot = xi_1_dot;
xi = [x_dot;W;W_R;W_L;W_A;];
I_g = eye(15, 15);
[K_tilde_g] = KK_tilde(INSECT, R, Q_R, Q_L, Q_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot);
[M_g, M_xi] = M(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
[M_tilde_g, M_tilde_xi] = M_tilde(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
L_tilde_g = M_g - 0.5*M_tilde_g;
L_tilde_xi = M_xi - 0.5*M_tilde_xi;
J_g_xi = JJ*xi;
A_tilde_g = blkdiag(zeros(3,3), hat(J_g_xi(4:6)), hat(J_g_xi(7:9)), hat(J_g_xi(10:12)), hat(J_g_xi(13:15)));
[F_g] = F_all(INSECT,x,R,Q_R,Q_L,Q_A,F_R,F_L,zeros(3,1),tau(4:6),tau(7:9),tau(10:12));

F_linear(1:15,1:15) = co_ad;
F_linear(1:15,16:30) = I_g;
F_linear(16:30,1:15) = JJ \ (-K_tilde_g + co_ad*KK - L_tilde_g + F_g);
F_linear(16:30,16:30) = JJ \ (A_tilde_g + co_ad*JJ - L_tilde_xi - LL);

% delta(:, k+1) = delta(:, k) + dt*F_linear*delta(:, k);
% trace_integ = trace_integ + trace(F_linear) * dt;

% xi=[xi_1;xi_2];
% xi_dot=JJ\( co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau);
% disp(norm(xi_dot - [xi_1_dot; xi_2_dot]));

% X_dot=[xi_1; xi_1_dot; F_linear*delta];
X_dot=[xi_1; xi_1_dot; reshape(F_linear*delta_mat, 900, 1);];
end

%% Functions

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

function [KK_til] = KK_tilde(INSECT, R, Q_R, Q_L, Q_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot)
KK_til_R = KK_tilde_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R, Q_R, x_ddot, W_dot, W_R_dot);
KK_til_L = KK_tilde_sub(INSECT.m_L, INSECT.mu_L, INSECT.nu_L, INSECT.J_L, R, Q_L, x_ddot, W_dot, W_L_dot);
KK_til_A = KK_tilde_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, x_ddot, W_dot, W_A_dot);

KK_til=zeros(15,15);
KK_til(1:3,4:6) = KK_til_R(1:3,4:6) + KK_til_L(1:3,4:6) + KK_til_A(1:3,4:6);
KK_til(1:3,7:9) = KK_til_R(1:3,7:9);
KK_til(1:3,10:12) = KK_til_L(1:3,7:9);
KK_til(1:3,13:15) = KK_til_A(1:3,7:9);

KK_til(4:6,4:6) = KK_til_R(4:6,4:6) + KK_til_L(4:6,4:6) + KK_til_A(4:6,4:6);
KK_til(4:6,7:9) = KK_til_R(4:6,7:9);
KK_til(4:6,10:12) = KK_til_L(4:6,7:9);
KK_til(4:6,13:15) = KK_til_A(4:6,7:9);

KK_til(7:9,4:6) = KK_til_R(7:9,4:6);
KK_til(7:9,7:9) = KK_til_R(7:9,7:9);

KK_til(10:12,4:6) = KK_til_L(7:9,4:6);
KK_til(10:12,10:12) = KK_til_L(7:9,7:9);

KK_til(13:15,4:6) = KK_til_A(7:9,4:6);
KK_til(13:15,13:15) = KK_til_A(7:9,7:9);
end

function [KK_til] = KK_tilde_sub(m, mu, nu, J, R, Q, x_ddot, W_dot, W_i_dot)
KK_til=zeros(9,9);
KK_til(1:3,4:6) = m*R*hat((hat(mu)+hat(Q*nu))*W_dot) + m*R*hat(Q*hat(nu)*W_i_dot);
KK_til(1:3,7:9) = -m*R*hat(W_dot)*Q*hat(nu) + m*R*Q*hat(hat(nu)*W_i_dot);
KK_til(4:6,4:6) = m*(hat(mu)+hat(Q*nu))*hat(R'*x_ddot);
KK_til(4:6,7:9) = m*hat(R'*x_ddot)*Q*hat(nu) - Q*hat(J*Q'*W_dot) + Q*J*hat(Q'*W_dot) ...
    -m*hat(mu)*hat(W_dot)*Q*hat(nu) - m* hat(hat(mu)*W_dot)*Q*hat(nu) ...
    -Q*hat(J*W_i_dot) + m*hat(mu)*Q*hat(hat(nu)*W_i_dot);
KK_til(7:9,4:6) = m*hat(nu)*Q'*hat(R'*x_ddot);
KK_til(7:9,7:9) = m*hat(nu)*hat(Q'*R'*x_ddot) + J*hat(Q'*W_dot) - m*hat(nu)*hat(Q'*hat(mu)*W_dot);
end

function [M_g, M_xi] = M(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
[M_g_R, M_xi_R] = M_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R, Q_R, x_dot, W, W_R);
[M_g_L, M_xi_L] = M_sub(INSECT.m_L, INSECT.mu_L, INSECT.nu_L, INSECT.J_L, R, Q_L, x_dot, W, W_L);
[M_g_A, M_xi_A] = M_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, x_dot, W, W_A);

M_g=zeros(15,15);
M_g(1:3,4:6) = M_g_R(1:3,4:6) + M_g_L(1:3,4:6) + M_g_A(1:3,4:6);
M_g(1:3,7:9) = M_g_R(1:3,7:9);
M_g(1:3,10:12) = M_g_L(1:3,7:9);
M_g(1:3,13:15) = M_g_A(1:3,7:9);

M_g(4:6,4:6) = M_g_R(4:6,4:6) + M_g_L(4:6,4:6) + M_g_A(4:6,4:6);
M_g(4:6,7:9) = M_g_R(4:6,7:9);
M_g(4:6,10:12) = M_g_L(4:6,7:9);
M_g(4:6,13:15) = M_g_A(4:6,7:9);

M_g(7:9,4:6) = M_g_R(7:9,4:6);
M_g(7:9,7:9) = M_g_R(7:9,7:9);

M_g(10:12,4:6) = M_g_L(7:9,4:6);
M_g(10:12,10:12) = M_g_L(7:9,7:9);

M_g(13:15,4:6) = M_g_A(7:9,4:6);
M_g(13:15,13:15) = M_g_A(7:9,7:9);

M_xi=zeros(15,15);
M_xi(1:3,1:3) = M_xi_R(1:3,1:3) + M_xi_L(1:3,1:3) + M_xi_A(1:3,1:3);
M_xi(1:3,4:6) = M_xi_R(1:3,4:6) + M_xi_L(1:3,4:6) + M_xi_A(1:3,4:6);
M_xi(1:3,7:9) = M_xi_R(1:3,7:9);
M_xi(1:3,10:12) = M_xi_L(1:3,7:9);
M_xi(1:3,13:15) = M_xi_A(1:3,7:9);

M_xi(4:6,1:3) = M_xi(1:3,4:6)';
M_xi(4:6,4:6) = M_xi_R(4:6,4:6) + M_xi_L(4:6,4:6) + M_xi_A(4:6,4:6);
M_xi(4:6,7:9) = M_xi_R(4:6,7:9);
M_xi(4:6,10:12) = M_xi_L(4:6,7:9);
M_xi(4:6,13:15) = M_xi_A(4:6,7:9);

M_xi(7:9,1:3) = M_xi(1:3,7:9)';
M_xi(7:9,4:6) = M_xi(4:6,7:9)';
M_xi(7:9,7:9) = M_xi_R(7:9,7:9);

M_xi(10:12,1:3) = M_xi(1:3,10:12)';
M_xi(10:12,4:6) = M_xi(4:6,10:12)';
M_xi(10:12,10:12) = M_xi_L(7:9,7:9);

M_xi(13:15,1:3) = M_xi(1:3,13:15)';
M_xi(13:15,4:6) = M_xi(4:6,13:15)';
M_xi(13:15,13:15) = M_xi_A(7:9,7:9);
end

function [M_g, M_xi] = M_sub(m, mu, nu, J, R, Q, x_dot, W, W_i)
M_g=zeros(9,9);
M_xi=zeros(9,9);

M_g(1:3,4:6) = -m*R*hat(hat((hat(mu)+hat(Q*nu))*W + (Q*hat(nu)*W_i))*W + (-hat(W)*Q*hat(nu)+Q*hat(hat(nu)*W_i))*W_i);
M_g(1:3,7:9) = -m*R*(hat(W)*(hat(W)*Q*hat(nu)-Q*hat(hat(nu)*W_i)) - hat(W)*Q*hat(hat(nu)*W_i) + Q*hat(hat(hat(nu)*W_i)*W_i));
M_g(4:6,4:6) = -m*((hat(mu)+hat(Q*nu))*hat(W)*hat(R'*x_dot) + hat(Q*hat(nu)*W_i)*hat(R'*x_dot));
M_g(4:6,7:9) = m*hat(hat(R'*x_dot)*W)*Q*hat(nu) - m*hat(R'*x_dot)*Q*hat(hat(nu)*W_i) + Q*hat(hat(J*Q'*W)*W_i) + ...
    Q*hat(W_i)*J*hat(Q'*W) - Q*hat(J*hat(Q'*W)*W_i) -Q*J*hat(W_i)*hat(Q'*W) + m*hat(mu)*hat(W)*Q*hat(hat(nu)*W_i) + ...
    m* hat(hat(mu)*W)*Q*hat(hat(nu)*W_i) + Q*hat(hat(J*W_i)*W_i) - m*hat(mu)*Q*hat(hat(hat(nu)*W_i)*W_i);
M_g(7:9,4:6) = -m*hat(nu)*Q'*hat(W)*hat(R'*x_dot) - m*hat(nu)*hat(W_i)*Q'*hat(R'*x_dot);
M_g(7:9,7:9) = m*hat(nu)*hat(Q'*hat(R'*x_dot)*W) - m*hat(nu)*hat(W_i)*hat(Q'*R'*x_dot) - J*hat(W_i)*hat(Q'*W) + m*hat(nu)*hat(W_i)*hat(Q'*hat(mu)*W);

M_xi(1:3,4:6) = m*R*(-hat(W)*(hat(mu)+hat(Q*nu)) + hat(Q*hat(nu)*W_i));
M_xi(1:3,7:9) = -m*R*(hat(W)*Q*hat(nu) + Q*hat(W_i)*hat(nu));
M_xi(4:6,4:6) = Q*hat(W_i)*(J*Q') - Q*J*hat(W_i)*Q' + m*hat(mu)*hat(Q*hat(nu)*W_i) + m* hat(Q*hat(nu)*W_i)*hat(mu);
M_xi(4:6,7:9) = Q*hat(W_i)*J - m*hat(mu)*Q*hat(W_i)*hat(nu);
M_xi(4:6,1:3) = M_xi(1:3,4:6)';
M_xi(7:9,1:3) = M_xi(1:3,7:9)';
M_xi(7:9,4:6) = M_xi(4:6,7:9)';
end

function [M_tilde_g, M_tilde_xi] = M_tilde(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
[M_tilde_g_R, M_tilde_xi_R] = M_tilde_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R, Q_R, x_dot, W, W_R);
[M_tilde_g_L, M_tilde_xi_L] = M_tilde_sub(INSECT.m_L, INSECT.mu_L, INSECT.nu_L, INSECT.J_L, R, Q_L, x_dot, W, W_L);
[M_tilde_g_A, M_tilde_xi_A] = M_tilde_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, x_dot, W, W_A);

M_tilde_g=zeros(15,15);

M_tilde_g(4:6,4:6) = M_tilde_g_R(4:6,4:6) + M_tilde_g_L(4:6,4:6) + M_tilde_g_A(4:6,4:6);
M_tilde_g(4:6,7:9) = M_tilde_g_R(4:6,7:9);
M_tilde_g(4:6,10:12) = M_tilde_g_L(4:6,7:9);
M_tilde_g(4:6,13:15) = M_tilde_g_A(4:6,7:9);

M_tilde_g(7:9,4:6) = M_tilde_g(4:6,7:9)';
M_tilde_g(7:9,7:9) = M_tilde_g_R(7:9,7:9);

M_tilde_g(10:12,4:6) = M_tilde_g(4:6,10:12)';
M_tilde_g(10:12,10:12) = M_tilde_g_L(7:9,7:9);

M_tilde_g(13:15,4:6) = M_tilde_g(4:6,13:15)';
M_tilde_g(13:15,13:15) = M_tilde_g_A(7:9,7:9);

M_tilde_xi=zeros(15,15);

M_tilde_xi(4:6,1:3) = M_tilde_xi_R(4:6,1:3) + M_tilde_xi_L(4:6,1:3) + M_tilde_xi_A(4:6,1:3);
M_tilde_xi(4:6,4:6) = M_tilde_xi_R(4:6,4:6) + M_tilde_xi_L(4:6,4:6) + M_tilde_xi_A(4:6,4:6);
M_tilde_xi(4:6,7:9) = M_tilde_xi_R(4:6,7:9);
M_tilde_xi(4:6,10:12) = M_tilde_xi_L(4:6,7:9);
M_tilde_xi(4:6,13:15) = M_tilde_xi_A(4:6,7:9);

M_tilde_xi(7:9,1:3) = M_tilde_xi_R(7:9,1:3);
M_tilde_xi(7:9,4:6) = M_tilde_xi_R(7:9,4:6);
M_tilde_xi(7:9,7:9) = M_tilde_xi_R(7:9,7:9);

M_tilde_xi(10:12,1:3) = M_tilde_xi_L(7:9,1:3);
M_tilde_xi(10:12,4:6) = M_tilde_xi_L(7:9,4:6);
M_tilde_xi(10:12,10:12) = M_tilde_xi_L(7:9,7:9);

M_tilde_xi(13:15,1:3) = M_tilde_xi_A(7:9,1:3);
M_tilde_xi(13:15,4:6) = M_tilde_xi_A(7:9,4:6);
M_tilde_xi(13:15,13:15) = M_tilde_xi_A(7:9,7:9);

end

function [M_tilde_g, M_tilde_xi] = M_tilde_sub(m, mu, nu, J, R, Q, x_dot, W, W_i)
M_tilde_g=zeros(9,9);
M_tilde_xi=zeros(9,9);

M_tilde_g(4:6,4:6) = -2*m*hat((hat(mu)+hat(Q*nu))*W + (Q*hat(nu)*W_i))*hat(R'*x_dot);
M_tilde_g(4:6,7:9) = 2*m*hat(R'*x_dot)*(hat(W)*Q*hat(nu)-Q*hat(hat(nu)*W_i))    ;
M_tilde_g(7:9,4:6) = M_tilde_g(4:6,7:9)';
M_tilde_g(7:9,7:9) = 2*( -m*(hat(nu)*hat(Q'*hat(W)*R'*x_dot) + hat(hat(nu)*W_i)*hat(Q'*R'*x_dot)) +...
     hat(J*Q'*W)*hat(Q'*W) - hat(Q'*W)*J*hat(Q'*W) + m*hat(nu)*hat(Q'*hat(W)*hat(mu)*W) + ...
     hat(J*W_i)*hat(Q'*W) + m*hat(hat(nu)*W_i)*hat(Q'*hat(mu)*W) );

M_tilde_xi(4:6,1:3) = -m*hat((hat(mu)+hat(Q*nu))*W)*R' - m*hat(Q*hat(nu)*W_i)*R';
M_tilde_xi(4:6,4:6) = m*hat(R'*x_dot)*(hat(mu)+hat(Q*nu));
M_tilde_xi(4:6,7:9) = m*hat(R'*x_dot)*Q*hat(nu);
M_tilde_xi(7:9,1:3) = -m*hat(nu)*Q'*hat(W)*R' - m*hat(hat(nu)*W_i)*Q'*R';
M_tilde_xi(7:9,4:6) = m*hat(nu)*Q'*hat(R'*x_dot) - hat(Q'*W)*(J*Q') + hat(J*Q'*W)*Q' - ...
    m*hat(nu)*Q'*hat(hat(mu)*W) +m*hat(nu)*Q'*hat(W)*hat(mu) +hat(J*W_i)*Q' +m*hat(hat(nu)*W_i)*Q'*hat(mu);
M_tilde_xi(7:9,7:9) = m*hat(Q'*R'*x_dot)*hat(nu) -hat(Q'*W)*J -m*hat(Q'*hat(mu)*W)*hat(nu);
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

function [F_g] = F_all(INSECT,x,R,Q_R,Q_L,Q_A,F_R,F_L,F_A,tau_R,tau_L,tau_A)
e3=[0 0 1]';

mg_R=INSECT.m_R*INSECT.g;
mg_L=INSECT.m_L*INSECT.g;
mg_A=INSECT.m_A*INSECT.g;
nu_R = INSECT.nu_R;
nu_L = INSECT.nu_L;
nu_A = INSECT.nu_A;
mu_R = INSECT.mu_R;
mu_L = INSECT.mu_L;
mu_A = INSECT.mu_A;

F_grav = zeros(15,15);
F_grav(4:6,4:6) = mg_R*hat(mu_R + Q_R*nu_R)*hat(R'*e3) + mg_L*hat(mu_L + Q_L*nu_L)*hat(R'*e3) ...
    + mg_A*hat(mu_A + Q_A*nu_A)*hat(R'*e3);
F_grav(4:6,7:9) = mg_R*hat(R'*e3)*Q_R*hat(nu_R);
F_grav(4:6,10:12) = mg_L*hat(R'*e3)*Q_L*hat(nu_L);
F_grav(4:6,13:15) = mg_A*hat(R'*e3)*Q_A*hat(nu_A);

F_grav(7:9,4:6) = F_grav(4:6,7:9)';
F_grav(7:9,7:9) = mg_R*hat(nu_R)*hat(Q_R'*R'*e3);

F_grav(10:12,4:6) = F_grav(4:6,10:12)';
F_grav(10:12,10:12) = mg_L*hat(nu_L)*hat(Q_L'*R'*e3);

F_grav(13:15,4:6) = F_grav(4:6,13:15)';
F_grav(13:15,13:15) = mg_A*hat(nu_A)*hat(Q_A'*R'*e3);

F_a = zeros(15,15);
F_a(1:3,4:6) = -R*hat(Q_R*F_R) -R*hat(Q_L*F_L) -R*hat(Q_A*F_A);
F_a(1:3,7:9) = -R*Q_R*hat(F_R);
F_a(1:3,10:12) = -R*Q_L*hat(F_L);
F_a(1:3,13:15) = -R*Q_A*hat(F_A);

F_a(4:6,7:9) = -hat(mu_R)*Q_R*hat(F_R);
F_a(4:6,10:12) = -hat(mu_L)*Q_L*hat(F_L);
F_a(4:6,13:15) = -hat(mu_A)*Q_A*hat(F_A);

F_tau = zeros(15,15);
F_tau(7:9,7:9) = hat(Q_R'*tau_R);
F_tau(10:12,10:12) = hat(Q_L'*tau_L);
F_tau(13:15,13:15) = hat(Q_A'*tau_A);

F_g = F_grav + F_a + F_tau;
end