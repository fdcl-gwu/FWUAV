%% Tests to verify the analyitical expressions of linearized dynamics
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');

load('sim_QS_x_hover.mat');
etas = [[1;0;0], [0;1;0], [0;0;1], [2;-1;1]];
epsilons = [1e-6, 1e-8];

for i=1:size(etas,2)
    for j=1:size(epsilons)
        eta = etas(:,i) * epsilons(j);
%         [lhs, rhs, output] = test_M_sub(epsilons(j), eta, N, INSECT, R, Q_R, x_dot, W, W_R);
%         [lhs, rhs, output] = test_M_tilde(epsilons(j), eta, N, INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
%         [lhs, rhs, output] = test_KK_tilde(epsilons(j), eta, N, INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot);
        [lhs, rhs, output] = test_full_eqn(epsilons(j), eta, N, INSECT, x, R, Q_R, Q_L, Q_A,...
            x_dot, W, W_R, W_L, W_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot, F_R, F_L, M_R, M_L, tau);
    end
end

function [lhs, rhs, output] = test_full_eqn(epsilon, eta, N, INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot, F_R, F_L, M_R, M_L, tau)
%%
for k=1:N
    % Verification of full linearization
    [JJ, EL_rhs, KK, LL, co_ad] = EL_equation_terms(INSECT, x(:, k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k), W(:, k), W_R(:, k), W_L(:, k), W_A(:, k), x_ddot(:, k), W_dot(:, k), W_R_dot(:, k), W_L_dot(:, k), W_A_dot(:, k), F_R(:, k), F_L(:, k), M_R(:, k), M_L(:, k), tau(:, k));
    [JJ_new, EL_rhs_new] = EL_equation_terms(INSECT, x(:, k), R(:,:,k)*expmso3(eta), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k), W(:, k), W_R(:, k), W_L(:, k), W_A(:, k), x_ddot(:, k), W_dot(:, k), W_R_dot(:, k), W_L_dot(:, k), W_A_dot(:, k), F_R(:, k), F_L(:, k), M_R(:, k), M_L(:, k), tau(:, k));
    delta_g = [zeros(3,1); eta; zeros(3,1); zeros(3,1); zeros(3,1);];
    delta_xi = [zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1);];
    
    xi = [x_dot(:, k); W(:, k); W_R(:, k); W_L(:, k); W_A(:, k);];
    xi_dot = [x_ddot(:, k); W_dot(:, k); W_R_dot(:, k); W_L_dot(:, k); W_A_dot(:, k);];
    [K_tilde_g] = KK_tilde(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_ddot(:, k), W_dot(:, k), W_R_dot(:, k), W_L_dot(:, k), W_A_dot(:, k));
    [M_g, M_xi] = M(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k), W(:, k), W_R(:, k), W_L(:, k), W_A(:, k));
    [M_tilde_g, M_tilde_xi] = M_tilde(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k), W(:, k), W_R(:, k), W_L(:, k), W_A(:, k));
    L_tilde_g = M_g - 0.5*M_tilde_g;
    L_tilde_xi = M_xi - 0.5*M_tilde_xi;
    J_g_xi = JJ*xi;
    A_tilde_g = blkdiag(zeros(3,3), hat(J_g_xi(4:6)), hat(J_g_xi(7:9)), hat(J_g_xi(10:12)), hat(J_g_xi(13:15)));
    [F_g] = F_all(INSECT,R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),F_R(:, k),F_L(:, k),zeros(3,1),tau(4:6,k),tau(7:9,k),tau(10:12,k));
    
    lhs = (EL_rhs_new - JJ_new*xi_dot) - (EL_rhs - JJ*xi_dot);
    rhs = ((-K_tilde_g + co_ad*KK - L_tilde_g + F_g)) * delta_g + ((A_tilde_g + co_ad*JJ - L_tilde_xi - LL)) * delta_xi;
    
    % Verification of linearization of position subsection
%     [JJ, EL_rhs, KK, LL, co_ad] = EL_equation_terms(INSECT, x(:, k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k), W(:, k), W_R(:, k), W_L(:, k), W_A(:, k), x_ddot(:, k), W_dot(:, k), W_R_dot(:, k), W_L_dot(:, k), W_A_dot(:, k), F_R(:, k), F_L(:, k), M_R(:, k), M_L(:, k), tau(:, k));
%     [JJ_new, EL_rhs_new] = EL_equation_terms(INSECT, x(:, k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:, k)+eta, W(:, k), W_R(:, k), W_L(:, k), W_A(:, k), x_ddot(:, k), W_dot(:, k), W_R_dot(:, k), W_L_dot(:, k), W_A_dot(:, k), F_R(:, k), F_L(:, k), M_R(:, k), M_L(:, k), tau(:, k));
%     lhs = EL_rhs_new(1:3) - EL_rhs(1:3);
%     [d_L_R d_L_L d_D_R d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R(:, k), W_L(:, k), W_R_dot(:, k), W_L_dot(:, k), x_dot(:, k), R(:,:,k), W(:, k), Q_R(:,:,k), Q_L(:,:,k), eta);
%     d_F_R=d_L_R+d_D_R;
%     d_F_L=d_L_L+d_D_L;
%     rhs = R(:,:,k) * (Q_R(:,:,k)*d_F_R + Q_L(:,:,k)*d_F_L);
    
    output = (lhs - rhs) ./ max(abs(lhs),abs(rhs));
    e = max(eps, epsilon^2);
    output((isnan(output) & rhs == 0) | (abs(lhs) < e & abs(rhs) < e)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [JJ, EL_rhs, KK, LL, co_ad] = EL_equation_terms(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot, F_R, F_L, M_R, M_L, tau)
%%
% This is the right method since even the aerodynamic forces depend on \xi.
% Do this to get more accurate linearization
% [L_R, L_L, D_R, D_L, M_R, M_L, ...
%     F_rot_R, F_rot_L, M_rot_R, M_rot_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
% F_R=L_R+D_R+F_rot_R;
% F_L=L_L+D_L+F_rot_L;
% M_R=M_R+M_rot_R;
% M_L=M_L+M_rot_L;
M_A=zeros(3,1);

f_a=[R*Q_R*F_R + R*Q_L*F_L; hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R; M_L; M_A];
f_a_1=f_a(1:3);
f_a_2=f_a(4:15);

% gravitational force and moment
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:3);
f_g_2=f_g(4:15);

xi_1=[x_dot]; 
xi_2=[W; W_R; W_L; W_A];
xi_2_dot=[W_dot; W_R_dot; W_L_dot; W_A_dot];

[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

% [JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_3_12(JJ);
% [LL_11, LL_12, LL_21, LL_22] = inertia_sub_decompose_3_12(LL);
% [~, ~, ~, co_ad_22] = inertia_sub_decompose_3_12(co_ad);
% 
% xi_1_dot = JJ_11\( -JJ_12*xi_2_dot -LL_11*xi_1 - LL_12*xi_2 + f_a_1 + f_g_1);
% f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1 + JJ_22*xi_2) ...
%     + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
f_tau = [zeros(3,1); -tau(4:6)-tau(7:9)-tau(10:12); Q_R'*tau(4:6); Q_L'*tau(7:9); Q_A'*tau(10:12);];

xi = [xi_1; xi_2;];
EL_rhs = co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau;

end

function [lhs, rhs, output] = test_KK_tilde(epsilon, eta, N, INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot)
%%
for k=1:N
    KK_til = KK_tilde(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_ddot(:,k), W_dot(:,k), W_R_dot(:,k), W_L_dot(:,k), W_A_dot(:,k));
    [J, ~] = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    [J_new, ~] = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k)*expmso3(eta), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    del_J = J_new - J;
    lhs = del_J* [x_ddot(:,k);W_dot(:,k);W_R_dot(:,k);W_L_dot(:,k);W_A_dot(:,k);];
    rhs = KK_til*[zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1); eta;];
    output = (lhs - rhs) ./ max(abs(lhs),abs(rhs));
    e = max(eps, epsilon^2);
    output((isnan(output) & rhs == 0) | (abs(lhs) < e & abs(rhs) < e)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [lhs, rhs, output] = test_M_tilde(epsilon, eta, N, INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
%%
for k=1:N
    [M_tilde_g, M_tilde_xi] = M_tilde(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    [~, K] = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    [~, K_new] = inertia(INSECT, R(:,:,k)*expmso3(eta), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    del_K = K_new - K;
    lhs = del_K'* [x_dot(:,k);W(:,k);W_R(:,k);W_L(:,k);W_A(:,k);];
    rhs = M_tilde_g*[zeros(3,1); eta; zeros(3,1); zeros(3,1); zeros(3,1);] + M_tilde_xi*[zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);];
    output = (lhs - rhs) ./ max(abs(lhs),abs(rhs));
    e = max(eps, epsilon^2);
    output((isnan(output) & rhs == 0) | (abs(lhs) < e & abs(rhs) < e)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [lhs, rhs, output] = test_M(epsilon, eta, N, INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
%%
for k=1:N
    [M_g, M_xi] = M(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    [~, K] = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    [~, K_new] = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k)+eta);
    del_K = K_new - K;
    lhs = del_K* [x_dot(:,k);W(:,k);W_R(:,k);W_L(:,k);W_A(:,k);];
    rhs = M_g*[zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1); zeros(3,1);] + M_xi*[zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);eta;];
    output = (lhs - rhs) ./ max(abs(lhs),abs(rhs));
    e = max(eps, epsilon^2);
    output((isnan(output) & rhs == 0) | (abs(lhs) < e & abs(rhs) < e)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [lhs, rhs, output] = test_M_tilde_sub(epsilon, eta, N, INSECT, R, Q_R, x_dot, W, W_R)
%%
for k=1:N
    [M_tilde_g, M_tilde_xi] = M_tilde_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k), x_dot(:,k), W(:,k), W_R(:,k));
    [~, K_i] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k), x_dot(:,k), W(:,k), W_R(:,k));
    [~, K_i_new] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k), x_dot(:,k), W(:,k)+eta, W_R(:,k));
    del_K_i = K_i_new - K_i;
    lhs = del_K_i'* [x_dot(:,k);W(:,k);W_R(:,k);];
    rhs = M_tilde_g*[zeros(3,1); zeros(3,1); zeros(3,1);] + M_tilde_xi*[zeros(3,1);eta;zeros(3,1);];
    output = (lhs - rhs) ./ abs(lhs);
    output((isnan(output) & rhs == 0) | (abs(lhs) < eps & abs(rhs) < eps)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [lhs, rhs, output] = test_M_sub(epsilon, eta, N, INSECT, R, Q_R, x_dot, W, W_R)
%%
for k=1:N
    [M_g, M_xi] = M_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k), x_dot(:,k), W(:,k), W_R(:,k));
    [~, K_i] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k), x_dot(:,k), W(:,k), W_R(:,k));
    [~, K_i_new] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.nu_R, INSECT.J_R, R(:,:,k), Q_R(:,:,k)*expmso3(eta), x_dot(:,k), W(:,k), W_R(:,k));
    del_K_i = K_i_new - K_i;
    lhs = del_K_i* [x_dot(:,k);W(:,k);W_R(:,k);];
    rhs = M_g*[zeros(3,1); zeros(3,1); eta;] + M_xi*[zeros(3,1);zeros(3,1);zeros(3,1);];
    output = (lhs - rhs) ./ abs(lhs);
    output((isnan(output) & rhs == 0) | (abs(lhs) < eps & abs(rhs) < eps)) = 0;
    if ~all(output < 1e-2)
        fprintf('There is an error\n');
        disp(lhs)
        disp(output)
    end
end
end

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom(INSECT, WK_R, WK_L, t, X)
%%
x=X(1:3);
x_dot=X(4:6);
% delta=X(7:12);
delta_mat=reshape(X(7:42), 6, 6);

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

% F_linear = zeros(30, 30);
% % trace_integ = 0;
% x_ddot = xi_1_dot;
% xi = [x_dot;W;W_R;W_L;W_A;];
% I_g = eye(15, 15);
% [K_tilde_g] = KK_tilde(INSECT, R, Q_R, Q_L, Q_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot);
% [M_g, M_xi] = M(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
% [M_tilde_g, M_tilde_xi] = M_tilde(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
% L_tilde_g = M_g - 0.5*M_tilde_g;
% L_tilde_xi = M_xi - 0.5*M_tilde_xi;
% J_g_xi = JJ*xi;
% A_tilde_g = blkdiag(zeros(3,3), hat(J_g_xi(4:6)), hat(J_g_xi(7:9)), hat(J_g_xi(10:12)), hat(J_g_xi(13:15)));
% [F_g] = F_all(INSECT,R,Q_R,Q_L,Q_A,F_R,F_L,zeros(3,1),tau(4:6),tau(7:9),tau(10:12));
% 
% F_linear(1:15,1:15) = co_ad;
% F_linear(1:15,16:30) = I_g;
% F_linear(16:30,1:15) = JJ \ (-K_tilde_g + co_ad*KK - L_tilde_g + F_g);
% F_linear(16:30,16:30) = JJ \ (A_tilde_g + co_ad*JJ - L_tilde_xi - LL);
% 
% X_dot=[xi_1; xi_1_dot; F_linear*delta];
% % X_dot=[xi_1; xi_1_dot; reshape(F_linear*delta_mat, 900, 1);];

% d_x_dot = delta(4:6);
% [d_L_R d_L_L d_D_R d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L, d_x_dot);
% d_F_R=d_L_R+d_D_R;
% d_F_L=d_L_L+d_D_L;
% delta_dot=zeros(6,1);
% delta_dot(1:3) = d_x_dot;
% delta_dot(4:6) = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L);
% X_dot=[xi_1; xi_1_dot; delta_dot];

F_linear = zeros(6, 6);
[d_L_R d_L_L d_D_R d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
d_F_R=d_L_R+d_D_R;
d_F_L=d_L_L+d_D_L;
F_linear(1:3, 4:6) = eye(3);
F_linear(4:6, 4:6) = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L);
X_dot=[xi_1; xi_1_dot; reshape(F_linear*delta_mat, 36, 1);];

end

function [KK_til] = KK_tilde(INSECT, R, Q_R, Q_L, Q_A, x_ddot, W_dot, W_R_dot, W_L_dot, W_A_dot)
%%
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
%%
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
%%
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
%%
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
%%
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
%%
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

function [F_g] = F_all(INSECT,R,Q_R,Q_L,Q_A,F_R,F_L,F_A,tau_R,tau_L,tau_A)
%%
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
