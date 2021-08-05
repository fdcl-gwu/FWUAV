function sim_QS_xR_hover
% simulate the position and attitude of thorax (x,R) for obtaining hover, 
% for given wing kinematics, abdomen attitude

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
load('morp_MONARCH', 'MONARCH');
INSECT=MONARCH;
filename='sim_QS_xR_hover';
load_past_data = true;

WK.f=10.2247;
% WK.beta=25.4292*pi/180;
WK.type='BermanWang';
WK.ab_type='varying';
WK.phi_max=75*pi/180;
WK.psi_max=5*pi/180;
WK.psi_N = 2; % or 1

N_single = 200;
N = 1 + 1*N_single;
x0=[0; 0; 0;];
final_pos = [0; 0; 0;];
% final_pos = [0.1; 0; -0.1/3;]; % Experimental trajectory

%% The optimization algorithm
if load_past_data
    load(filename, 'WK', 'WK_arr', 'solutions', 'output', 'lb', 'ub');
else
    A = []; b = []; Aeq = []; beq = [];
    % Initial value of WK_arr = [beta, phi_m, phi_K, phi_0, theta_m, theta_C, theta_0, theta_a, psi_m, psi_a, psi_0, x_dot1, x_dot2, x_dot3, theta_A_m, theta_A_0, theta_A_a, freq, theta_B02, W_B02]
    WK_arr0 = [-0.2400   0.7806  0.4012   0.7901    0.6981    2.9999    0.2680    0.1050    5*pi/180    0.9757    5*pi/180   -0.1000   -0.0000 -0.1000 0.0937 pi/15 0 WK.f 0 0];
	lb = [-pi/8, 0, 0, -pi/2, 0, 0, -pi/6, -pi/2, 0, -pi, -5*pi/180, -2, -2, -2, 0, pi/15, -pi, WK.f*(1-0.15), 0, -5];
	ub = [pi/5, pi/2, 1, pi/2, 40*pi/180, 3, pi/6, pi/2, 5*pi/180, pi, 5*pi/180, 2, 2, 2, pi/12, pi/4, pi, WK.f*(1+0.15), 45*pi/180, 5];
    % lb = [-pi/12, 0, 0, -pi/2, 0, 0, -pi/6, -pi/2, 0, -pi, -5*pi/180, -2, -2, -2, 0, pi/15, -pi, WK.f*(1-0.15), 0, -5];
    % ub = [pi/8, pi/2, 1, pi/2, 40*pi/180, 3, pi/6, pi/2, 5*pi/180, pi, 5*pi/180, 2, 2, 2, pi/12, pi/6, pi, WK.f*(1+0.15), 45*pi/180, 5];
    nonlcon = @(WK_arr) traj_condition(WK_arr, WK, INSECT, N, x0, final_pos);

    tic;
    rng default; % For reproducibility

    % % FMINCON
    % WK_arr0 = [0.7782    0.6355    0.2866   -0.6599    0.6893    2.1703    0.0098   -0.1410    0.0196    0.2506   -0.0003 -0.2458   -0.0000    0.0230    0.1970    0.4696    1.4270   11.6689  0.8319    1.4862];
    % options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
    %     'MaxFunctionEvaluations',5000,'PlotFcn',@optimplotfval,'UseParallel',true);
    % [WK_arr, fval, exitflag, output] = fmincon(@(WK_arr) objective_func(WK_arr, WK, INSECT, N, x0, final_pos),...
    %     WK_arr0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    % % MULTISTART, PARTICLESWARM
    ptmatrix(1, :) = [0.2950    0.8865    0.7242   -0.3732    0.5625    0.6328   -0.0129    0.3241    0.0419    0.2634   -0.0000    0.6994    0.0000   -0.6685 0.1648    0.2433    3.1416    9.3263 0 0];
    ptmatrix(2, :) = [-0.2799    0.7425    0.6669    0.5665    0.6981    3.0000    0.3914   -0.2456    0.0000   -2.8879    0.0524   -0.1000   -0.0000   -0.1000 0.2094   0.2669    1.5707   11.7584 0 0];
    ptmatrix(3, :) = [-0.3324    1.4397-0.6658    0.3889    0.6658    0.1137    1.7289   -0.5236    1.5708    5*pi/180   -1.0691    5*pi/180   -0.0441    0.0137 -0.0783 0 0.21 0 WK.f 0 0];
    ptmatrix(4, :) = [0.0577    pi/2-1.2217    0.3587    1.2217  0.5233    2.7760    0.2207    0.0059    5*pi/180    0.9429    0.0291 0 0 0 10*pi/180  0.21 0 WK.f 0 0];
    ptmatrix(5, :) = [0.6273    0.6355    0.2866   -0.6599    0.6893    2.1703    0.0098   -0.1410    0.0196    0.2506   -0.0003 -0.2458   -0.0000    0.0230    0.1970    0.4696    1.4270   11.6689  0.5319    1.4862];
    ptmatrix(6, :) = [0.3273    0.7399    0.9091    0.5642    0.0000    0.2846   -0.5236    1.1408    0.0244    2.2337    0.0278 0.2275    0.0000    0.1058    0.2093    0.2778    1.8785   10.9897   0.7069    1.8056];
    N_points = 50;
    ptmatrix(7:N_points, :) = lb + rand(N_points-6, length(WK_arr0)) .* (ub - lb);
    tpoints = CustomStartPointSet(ptmatrix);
    ms = MultiStart('Display','iter','PlotFcn',@gsplotbestf,'MaxTime',12*3600);
    options = optimoptions(@fmincon,'Algorithm','sqp',...
        'UseParallel',true,'MaxIterations',2000,'MaxFunctionEvaluations',6000); ...
        % 'ConstraintTolerance',1e-10,'StepTolerance',1e-8,'ObjectiveLimit',0);
    problem = createOptimProblem('fmincon','objective',@(WK_arr) objective_func(WK_arr, WK, INSECT, N, x0, final_pos),...
        'x0',WK_arr0,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
    [WK_arr, fval, exitflag, output, solutions] = run(ms, problem, tpoints);
    % WK_arr = [0.3273    0.7399    0.9091    0.5642    0.0000    0.2846   -0.5236    1.1408    0.0244    2.2337    0.0278 0.2275    0.0000    0.1058    0.2093    0.0778    1.8785   10.9897   -1.3069    1.8056];

    fprintf('Optimization has been completed\n');
    disp(output);
    toc;
end

sol_arr = [];
for i=1:length(solutions)
	if all(solutions(i).Output.bestfeasible.x == solutions(i).X) && (solutions(i).Output.lssteplength < 1)
		sol_arr = [sol_arr; i, solutions(i).Fval, solutions(i).Output.firstorderopt, ...
			solutions(i).Output.funcCount, solutions(i).Output.constrviolation];
	end
end
sol_arr = array2table(sol_arr, 'VariableNames', ...
    {'Index','Fval','Firstorderopt','FuncCount','Constrviolation'});
% sol_arr = sortrows(sol_arr, 'Firstorderopt');
disp(sol_arr);
if ~exist('sol_idx')
    sol_idx = sol_arr{1, 1};
end
% sol_idx = 1;
WK_arr = solutions(sol_idx).X;
fval = solutions(sol_idx).Fval;

%%
[WK, x_dot0, R0, W0] = get_WK(WK, WK_arr);
X0=[x0; reshape(R0,9,1); x_dot0; W0];
% load('sim_QS_xR_hover_control_opt_single.mat', 'X0'); N=201; T=2/WK.f;

N_periods=3;
N=1+N_periods*N_single; % minimum 30 / period
T=N_periods/WK.f;
t=linspace(0,T,N);

% [t, X]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
X = crgr_xR_mex(INSECT, WK, WK, t, X0);

x = X(:,1:3)';
x_dot = X(:,13:15)';

for k=1:N    
    [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),...
        theta_A(k), W(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
        W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
        M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k), ...
        Euler_R(:,k), Euler_R_dot(:,k)]= eom_QS_xR(INSECT, WK, WK, t(k), X(k,:)');
    F_B(:,k) = Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);
    theta_B(k) = R2axang(R(:,:,k));
end

x_ddot = X_dot(13:15,:);

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

end

function [c,ceq] = traj_condition(WK_arr, WK, INSECT, N, x0, final_pos)
%% Nonlinear constraint
[WK, x_dot0, R0, W0] = get_WK(WK, WK_arr);
X0=[x0; reshape(R0,9,1); x_dot0; W0];

T=1/WK.f;
t=linspace(0,T,N);

% [t, X]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
X = crgr_xR_mex(INSECT, WK, WK, t, X0);

c(1) = abs(WK.phi_0) + WK.phi_m - WK.phi_max;
c(2) = abs(WK.psi_0) + WK.psi_m - WK.psi_max;
% R=reshape(X(:,4:12)',3,3,N);
% for k=1:N
%     theta_B(k) = acos((trace(R(:,:,k))-1)/2);
% end
% theta_B_max = max(theta_B);
% c(3) = theta_B_max - WK.theta_B_max;
ceq(1:3) = X(1,1:3)' - (X(end,1:3)' -final_pos); % x
ceq(4:12) = 1e-2*(X(1,4:12)' - X(end,4:12)'); % R
ceq(13:15) = 1e-1*(X(1,13:15)' - X(end,13:15)'); % x_dot
ceq(16:18) = 1e-3*(X(1,16:18)' - X(end,16:18)'); % W
% ceq(1:3) = 10* (X(1,1:3)' - (X(end,1:3)' -final_pos)); % x
% ceq(4:12) = (X(1,4:12)' - X(end,4:12)'); % R
% ceq(13:15) = 5*(X(1,13:15)' - X(end,13:15)'); % x_dot
% ceq(16:18) = (X(1,16:18)' - X(end,16:18)'); % W
if any(isnan(ceq), 'all')
    ceq(1:18) = 1/eps;
end
end

function [J] = objective_func(WK_arr, WK, INSECT, N, x0, final_pos)
%% Objective function
[WK, x_dot0, R0, W0] = get_WK(WK, WK_arr);
X0=[x0; reshape(R0,9,1); x_dot0; W0];

T=1/WK.f;
t=linspace(0,T,N);

% [t, X]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
X = crgr_xR_mex(INSECT, WK, WK, t, X0);

x=X(:,1:3)';
x_dot=X(:,13:15)';
m=INSECT.m;
E = 0.5*m*(vecnorm(x_dot).^2) - m*9.81*x(3,:);
E_dot = diff(E')./diff(t');
E_dot = [E_dot; E_dot(end)];

gam = 1e3;
del = 1e4;
% J = gam * abs(trapz(t, E_dot)) + del * abs(trapz(t, E)); % new data
J = gam * trapz(t, abs(E_dot)) + del * trapz(t, abs(E));

% % Augmentation of constraints
% c(1) = abs(WK.phi_0) + WK.phi_m - WK.phi_max;
% c(2) = abs(WK.psi_0) + WK.psi_m - WK.psi_max;
% ceq(1:3) = X(1,1:3)' - (X(end,1:3)' -final_pos); % x
% ceq(4:12) = 1e-2*(X(1,4:12)' - X(end,4:12)'); % R
% ceq(13:15) = 1e-1*(X(1,13:15)' - X(end,13:15)'); % x_dot
% ceq(16:18) = 1e-3*(X(1,16:18)' - X(end,16:18)'); % W
% w_c = 1; w_ceq = 1e6;
% for k=1:length(c)
%     if c(k) >= 0
%         J = 1/eps;
%     else
%         J = J - w_c*log(-c(k));
%     end
% end
% J = J + w_ceq * (exp(norm(ceq)) - 1);

if isnan(J)
    J = 1/eps;
end

end

function [WK, x_dot0, R0, W0] = get_WK(WK, WK_arr)
WK.beta=WK_arr(1);
WK.phi_m=WK_arr(2);
WK.phi_K=WK_arr(3);
WK.phi_0=WK_arr(4);

WK.theta_m=WK_arr(5);
WK.theta_C=WK_arr(6);
WK.theta_0=WK_arr(7);
WK.theta_a=WK_arr(8);

WK.psi_m=WK_arr(9);
WK.psi_a=WK_arr(10);
WK.psi_0=WK_arr(11);

x_dot0=[WK_arr(12), WK_arr(13), WK_arr(14)]';

WK.theta_A_m = WK_arr(15);
WK.theta_A_0 = WK_arr(16);
WK.theta_A_a = WK_arr(17);
WK.f = WK_arr(18);

e2 = [0; 1; 0];
R0 = expmhat(WK_arr(19)*e2);
W0 = [0; WK_arr(20); 0];
end
