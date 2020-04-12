function sim_QS_x_hover
% simulate the position of thorax (x) for obtaining hover, 
% for given thorax attiude, wing kinematics, abdomen attitude

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
load('./sim_data/other_insects/morp_hawk', 'INSECT');
filename='sim_QS_x_hover_hawk';

WK.f=INSECT.f;
WK.type='BermanWang';
WK.ab_type='varying';
WK.bo_type='varying';
WK.phi_max=75*pi/180;
WK.psi_max=3*pi/180;
WK.psi_N = 2; % or 1

N=1001;
x0=[0 0 0]';
final_pos = [0; 0; 0;];
% final_pos = [0.1; 0; -0.1/3;]; % Experimental trajectory

%% The optimization algorithm
A = []; b = []; Aeq = []; beq = [];
% Initial value of WK_arr = [beta, phi_m, phi_K, phi_0, theta_m, theta_C, theta_0, theta_a, psi_m, psi_a, psi_0, x_dot1, x_dot2, x_dot3, theta_B_m, theta_B_0, theta_B_a, theta_A_m, theta_A_0, theta_A_a, freq]
WK_arr0 = [-0.2400   0.7806  0.4012   0.7901    0.6981    2.9999    0.2680    0.1050    8*pi/180    0.9757    5*pi/180   -0.1000   -0.0000 -0.1000 0 0 0 0.0937 0 0 WK.f];
lb = [-pi/2, 0, 0, -pi/2, 0, 0, -pi/6, -pi/2, 0, -pi, -5*pi/180, -1, -1, -1, 0, -pi/9, -pi, 0, -pi/4, -pi, WK.f*(1-0.15)];
ub = [pi/2, pi/2, 1, pi/2, 40*pi/180, 3, pi/6, pi/2, 5*pi/180, pi, 5*pi/180, 1, 1, 1, pi/18, pi/3, pi, pi/15, pi/4, pi, WK.f*(1+0.15)];
nonlcon = @(WK_arr) traj_condition(WK_arr, WK, INSECT, N, x0, final_pos);

tic;
rng default; % For reproducibility

% % FMINCON
% WK_arr0 = [ 0.2950    0.8865    0.7242   -0.3732    0.5625    0.6328   -0.0129    0.3241    0.0419    0.2634   -0.0000    0.6994    0.0000   -0.6685 0.1102    0.4683   -0.1631    0.1648    0.0433    3.1416    9.3263];
% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
%     'MaxFunctionEvaluations',5000,'PlotFcn',@optimplotfval,'UseParallel',true);
% [WK_arr, fval, exitflag, output] = fmincon(@(WK_arr) objective_func(WK_arr, WK, INSECT, N, x0, final_pos),...
%     WK_arr0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% % MULTISTART, PARTICLESWARM
ptmatrix(1, :) = [0.2950    0.8865    0.7242   -0.3732    0.5625    0.6328   -0.0129    0.3241    0.0419    0.2634   -0.0000    0.6994    0.0000   -0.6685 0.1102    0.4683   -0.1631    0.1648    0.0433    3.1416    9.3263];
ptmatrix(2, :) = [-0.2799    0.7425    0.6669    0.5665    0.6981    3.0000    0.3914   -0.2456    0.0000   -2.8879    0.0524   -0.1000   -0.0000   -0.1000 0.1745    1.0472    1.5708    0.2094   -0.1669    1.5707   11.7584];
ptmatrix(3, :) = [-0.4324    1.4397-0.6658    0.3889    0.6658    0.1137    1.7289   -0.5236    1.5708    8*pi/180   -1.0691    5*pi/180   -0.0441    0.0137 -0.0783 0 15*pi/180 0 10*pi/180  0 0 WK.f];
ptmatrix(4, :) = [0.0577    pi/2-1.2217    0.3587    1.2217  0.5233    2.7760    0.2207    0.0059    8*pi/180    0.9429    0.0291 0 0 0 0 15*pi/180 0 10*pi/180  0 0 WK.f];
ptmatrix(5, :) = [-1.1117    0.8052    0.9530    0.3491    0.0039    1.5236    0.5158    0.8012    0.1344   -1.5344    0.0771    0.0981    0.0000   -0.0814  0  0.3219 0 0.2869   -0.3374    1.5438   11.8923];
ptmatrix(6, :) = [0.1815    pi/2-1.2217   0.3218    -1.2217    0.5196    2.9411    0.1942    0.0157    8*pi/180    0.9508   -0.0505 0 0 0 0 15*pi/180 0 10*pi/180  0 0 WK.f];

ptmatrix(7:26, :) = lb + rand(20, length(WK_arr0)) .* (ub - lb);
tpoints = CustomStartPointSet(ptmatrix);
ms = MultiStart('Display','iter','PlotFcn',@gsplotbestf,'MaxTime',12*3600);
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'UseParallel',true);%'ConstraintTolerance',1e-5,'StepTolerance',1e-8,'OptimalityTolerance',1e-5);
problem = createOptimProblem('fmincon','objective',@(WK_arr) objective_func(WK_arr, WK, INSECT, N, x0, final_pos),...
    'x0',WK_arr0,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
[WK_arr, fval, exitflag, output, solutions] = run(ms, problem, tpoints);

% hybridoptions = optimoptions('fmincon','MaxFunctionEvaluations',5000,'UseParallel',true);
% options = optimoptions('particleswarm','PlotFcn',@pswplotbestf,...
%     'UseParallel',true,'InitialSwarmMatrix',ptmatrix,'MaxTime',12*3600);
% %     'HybridFcn',{@fmincon,hybridoptions});
% [WK_arr, fval, exitflag, output] = particleswarm(@(WK_arr) objective_func(WK_arr, WK, INSECT, N, x0, final_pos),...
%     length(WK_arr0), lb, ub, options);

fprintf('Optimization has been completed\n');
disp(output);
toc;

%%
[WK, x_dot0] = get_WK(WK, WK_arr);
X0=[x0; x_dot0];

N=1001;
T=3/WK.f;
t=linspace(0,T,N);

[t, X]=ode45(@(t,X) eom_QS_x(INSECT, WK, WK, t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x = X(:,1:3)';
x_dot = X(:,4:6)';

for k=1:N    
    [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), theta_B(k),...
        theta_A(k), W(:,k), W_dot(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
        W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
        M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k)]= eom_QS_x(INSECT, WK, WK, t(k), X(k,:)');
    F_B(:,k) = Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);    
    [Euler_R(:,k), Euler_R_dot(:,k), Euler_R_ddot(:,k)] = wing_kinematics(t(k),WK);
end

x_ddot = X_dot(4:6,:);

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
[WK, x_dot0] = get_WK(WK, WK_arr);
X0=[x0; x_dot0];

T=1/WK.f;
t=linspace(0,T,N);
[t, X]=ode45(@(t,X) eom_QS_x(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
c(1) = abs(WK.phi_0) + WK.phi_m - WK.phi_max;
c(2) = abs(WK.psi_0) + WK.psi_m - WK.psi_max;
ceq(1:3) = X(1,1:3)' - (X(end,1:3)' -final_pos);
ceq(4:6) = 0.1*(X(1,4:6)' - X(end,4:6)');
if any(isnan(ceq), 'all')
    ceq(1:6) = 1/eps;
end
end

function [J] = objective_func(WK_arr, WK, INSECT, N, x0, final_pos)
%% Objective function
[WK, x_dot0] = get_WK(WK, WK_arr);
X0=[x0; x_dot0];

T=1/WK.f;
t=linspace(0,T,N);
[t, X]=ode45(@(t,X) eom_QS_x(INSECT, WK, WK, t, X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,4:6)';
m=INSECT.m;
E = 0.5*m*(vecnorm(x_dot).^2) - m*9.81*x(3,:);
E_dot = diff(E')./diff(t);
E_dot = [E_dot; E_dot(end)];

gam = 1e3;
del = 1e4;
J = gam * trapz(t, abs(E_dot)) + del * trapz(t, abs(E));

% % Augmentation of constraints
% c(1) = abs(WK.phi_0) + WK.phi_m - WK.phi_max;
% c(2) = abs(WK.psi_0) + WK.psi_m - WK.psi_max;
% ceq(1:3) = x(:, 1) - (x(:, end) -final_pos);
% ceq(4:6) = 0.1*(x_dot(:, 1) - x_dot(:, end));
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

function [WK, x_dot0] = get_WK(WK, WK_arr)
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

WK.theta_B_m = WK_arr(15);
WK.theta_B_0 = WK_arr(16);
WK.theta_B_a = WK_arr(17);
WK.theta_A_m = WK_arr(18);
WK.theta_A_0 = WK_arr(19);
WK.theta_A_a = WK_arr(20);
WK.f = WK_arr(21);
end
