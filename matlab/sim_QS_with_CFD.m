function MAIN
clc
% simulate the position of thorax (x),
% for given thorax attiude, wing kinematics, abdomen attitude
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');

load('morp_MONARCH');
INSECT=MONARCH;

% Create wing kinematic information:
WK.f=10.2247;
WK.beta=25.4292*pi/180;
WK.type='Monarch';
WK.ab_type='experimental';
WK.bo_type='experimental';

% Set simulation parameters:
SIM.stepsPerPeriod = 100;
SIM.numPeriods = 10;
SIM.N = SIM.numPeriods*SIM.stepsPerPeriod+1;
SIM.T = SIM.numPeriods/WK.f;
t=linspace(0,SIM.T,SIM.N);
SIM.dt = t(2)-t(1);

% Intitial conditions:
x0=[0 0 0]';
x_dot0=[1.0 0 -0.3]';

% To use old CFD results, set both of these to 1.
% To run and use new CFD results, set cfdOnFlag = 1 and collectData = 0;
SIM.cfdOnFlag = 0;
SIM.collectData = 0;

SIM.dynamics = 'x'; % Choose either x or xR
% Note that xR runs into errors if run beyond 1 flapping period
SIM.QS_rot_force = 'off'; % Choose either rotational forces on or off here

if strcmp(SIM.dynamics,'x') == true
    if strcmp(SIM.QS_rot_force,'off') == true
        filename='sim_QS_x';
    elseif strcmp(SIM.QS_rot_force,'on') == true
        filename='sim_QSR_x';
    end
    X0=[x0; x_dot0];
    X = zeros(length(X0),length(t));
    XDot = zeros(length(X0),length(t));
    % Time integrate the states X:
    X = time_int_adam_bash(@eom_x,INSECT,WK,WK,t,X,XDot,X0,SIM);
    % Break out data for plotting later:
    % Note: x = X(:,1:3); x_dot = X(:,4:6);
    
    % Initialize final run data arrays:
    % 3x3xN arrays:
    R=zeros(3,3,SIM.N); Q_R=R; Q_L=R; Q_A=R;
    % 3xN arrays:
    W=zeros(3,SIM.N); W_R=W; W_L=W; W_A=W; F_R_CFD=W; F_I_CFD=W; F_B_CFD=W; 
    M_I_CFD=W; M_R_CFD=W; M_L_CFD=W; F_R_QS=W; F_L_QS=W; M_R_QS=W; M_L_QS=W;
    M_B_CFD=W; F_R=W; F_L=W; M_R=W; M_L=W; F_B_QS=W; F_I_QS=W;
    Euler_R=W; Euler_R_dot=W; Euler_R_ddot=W; 
    U_alpha_R_dot=W; U_alpha_L_dot=W; U_R=W; U_L=W;
    % 15xN arrays:
    f_a=zeros(15,SIM.N); f_g=f_a; f_tau=f_a;
    % 12xN arrays:
    tau=zeros(12,SIM.N);
    % 1xN arrays:
    theta_B=zeros(1,SIM.N); theta_A=theta_B;
    alpha_R=zeros(1,SIM.N); alpha_L=zeros(1,SIM.N);
    
    SIM.collectData = 1;
    for k=1:SIM.N
        % Time history info:
        [~,R(:,:,k),Q_R(:,:,k),Q_L(:,:,k),Q_A(:,:,k),theta_B(k), ...
            theta_A(k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k), ...
            F_I_CFD(:,k),F_B_CFD(:,k),F_R_CFD(:,k),...
            M_I_CFD(:,k),M_B_CFD(:,k),M_R_CFD(:,k),M_L_CFD(:,k),...
            F_R(:,k),F_L(:,k),M_R(:,k),M_L(:,k),f_a(:,k),f_g(:,k), ...
            f_tau(:,k),tau(:,k), ...
            alpha_R(:,k), alpha_L(:,k), ...
            U_alpha_R_dot(:,k), U_alpha_L_dot(:,k), U_R(:,k), U_L(:,k)] = ...
            eom_x(INSECT,WK,WK,t(k),X(:,k),k,SIM);
        
        F_B_QS(:,k)=Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);
        F_I_QS(:,k)=R(:,:,k)*F_B_QS(:,k);
        F_R_QS(:,k)=F_R(:,k);
        M_R_QS(:,k)=M_R(:,k);
        F_L_QS(:,k)=F_L(:,k);
        M_L_QS(:,k)=M_L(:,k);
        
        % Wing kinematics:
        [Euler_R(:,k),Euler_R_dot(:,k),Euler_R_ddot(:,k)] = ...
            wing_kinematics(t(k),WK);
    end
elseif strcmp(SIM.dynamics,'xR') == true
    if strcmp(SIM.QS_rot_force,'off') == true
        filename='sim_QS_xR';
    elseif strcmp(SIM.QS_rot_force,'on') == true
        filename='sim_QSR_xR';
    end
    % Additional initial conditions for unprescribed body dynamics
    theta_B0 = 20; % Initial body pitch
    R0=axang2rotm([0 1 0 deg2rad(theta_B0)]);
    W0=[0; 10; 0]; % Initial body rotation rates (x, y, and z)
    X0=[x0; reshape(R0,9,1); x_dot0; W0];
    
    X = zeros(length(X0),length(t));
    XDot = zeros(length(X0),length(t));
    % Time integrate the states X:
    X = time_int_adam_bash(@eom_xR,INSECT,WK,WK,t,X,XDot,X0,SIM);
    % NOTE: x=X(:,1:3)'; R=X(:,4:12)'; x_dot=X(:,13:15)'; W=X(:,16:18)';
    
    % Initialize final run data arrays:
    % 3x3xN arrays:
    R=zeros(3,3,SIM.N); Q_R=R; Q_L=R; Q_A=R;
    % 3xN arrays:
    W=zeros(3,SIM.N); W_R=W; W_L=W; W_A=W; F_R_CFD=W; F_I_CFD=W; F_B_CFD=W; 
    M_I_CFD=W; M_R_CFD=W; M_L_CFD=W; F_R_QS=W; F_L_QS=W; M_R_QS=W; M_L_QS=W;
    M_B_CFD=W; F_R=W; F_L=W; M_R=W; M_L=W; F_B_QS=W; F_I_QS=W;
    Euler_R=W; Euler_R_dot=W; Euler_R_ddot=W; 
    U_alpha_R_dot=W; U_alpha_L_dot=W; U_R=W; U_L=W;
    % 15xN arrays:
    f_a=zeros(15,SIM.N); f_g=f_a; f_tau=f_a;
    % 9xN arrays:
    tau=zeros(9,SIM.N);
    % 1xN arrays:
    theta_B=zeros(1,SIM.N); theta_A=theta_B;
    alpha_R=zeros(1,SIM.N); alpha_L=zeros(1,SIM.N);
    
    % This is always false at this step, since we just want to write out the
    % final arrays if we've used the CFD aero model:
    SIM.cfd_new_run_bool = false;
    
    for k=1:SIM.N
        % Time history info:
        [~,R(:,:,k),Q_R(:,:,k),Q_L(:,:,k),Q_A(:,:,k),theta_B(k), ...
            theta_A(k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k), ...
            F_I_CFD(:,k),F_B_CFD(:,k),F_R_CFD(:,k),...
            M_I_CFD(:,k),M_B_CFD(:,k),M_R_CFD(:,k),M_L_CFD(:,k),...
            F_R(:,k),F_L(:,k),M_R(:,k),M_L(:,k),f_a(:,k),f_g(:,k), ...
            f_tau(:,k),tau(:,k), ...
            alpha_R(:,k), alpha_L(:,k), ...
            U_alpha_R_dot(:,k), U_alpha_L_dot(:,k), U_R(:,k), U_L(:,k)] = ...
            eom_xR(INSECT,WK,WK,t(k),X(:,k),k,SIM);
        
        F_B_QS(:,k)=Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);
        F_I_QS(:,k)=R(:,:,k)*F_B_QS(:,k);
        F_R_QS(:,k)=F_R(:,k);
        M_R_QS(:,k)=M_R(:,k);
        F_L_QS(:,k)=F_L(:,k);
        M_L_QS(:,k)=M_L(:,k);
        
        % Wing kinematics:
        [Euler_R(:,k),Euler_R_dot(:,k),Euler_R_ddot(:,k)] = ...
            wing_kinematics(t(k),WK);
    end
end

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics 
% object
tosave = ...
    cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);
movefile([filename,'.mat'],'./sim_data')

% run('checkSol.m')
% run('POST/fig_comp_QSR_vs_CFD_x.m');
run('plotting/fig_comp_QSR_vs_CFD_x.m');

end
