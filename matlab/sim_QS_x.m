function sim_QS_x
% simulate the position of thorax (x), 
% for given thorax attiude, wing kinematics, abdomen attitude

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='sim_QS_x';

load('morp_MONARCH', 'MONARCH');
INSECT=MONARCH;

WK.f=10.2247;
WK.beta=25.4292*pi/180;
WK.type='Monarch';
WK.ab_type='experimental';
WK.bo_type='experimental';

N=5001;
T=5/WK.f;
t=linspace(0,T,N);

x0=[0 0 0]';
x_dot0=[1.0 0 -0.3]';

X0=[x0; x_dot0];

[t X]=ode45(@(t,X) eom_QS_x(INSECT, WK, WK, t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,4:6)';

R=zeros(3,3,N);
for k=1:N    
    [~, R(:,:,k) Q_R(:,:,k) Q_L(:,:,k) Q_A(:,:,k) theta_B(k) theta_A(k) W(:,k) W_dot(:,k) W_R(:,k) W_R_dot(:,k) W_L(:,k) W_L_dot(:,k) W_A(:,k) W_A_dot(:,k) F_R(:,k) F_L(:,k) M_R(:,k) M_L(:,k) f_a(:,k) f_g(:,k) f_tau(:,k) tau(:,k)]= eom_QS_x(INSECT, WK, WK, t(k), X(k,:)');
    F_B(:,k)=Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);    
    [Euler_R(:,k), Euler_R_dot(:,k), Euler_R_ddot(:,k)] = wing_kinematics(t(k),WK);
end

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
