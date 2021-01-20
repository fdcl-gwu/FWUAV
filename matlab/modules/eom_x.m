function [X_dot,R,Q_R,Q_L,Q_A,theta_B,theta_A,W,W_R,W_L,W_A,...
    F_I_CFD,F_B_CFD,F_R_CFD,M_I_CFD,M_B_CFD,M_R_CFD,M_L_CFD,...
    F_R,F_L,M_R,M_L,f_a,f_g,f_tau,tau,...
    alpha_R,alpha_L,U_alpha_R_dot,U_alpha_L_dot,U_R,U_L] = ...
    eom_x(INSECT,WK_R,WK_L,t,X,timeStep,SIM)

x=X(1:3);
x_dot=X(4:6);

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);

[Q_R,Q_L,W_R,W_L,W_R_dot,W_L_dot] = ...
    wing_attitude(WK_R.beta,Euler_R,Euler_L,Euler_R_dot,Euler_L_dot,...
    Euler_R_ddot,Euler_L_ddot);

[R,W,W_dot,theta_B] = body_attitude(t,WK_R.f,WK_R); %time-varying thorax
[Q_A,W_A,W_A_dot,theta_A] = abdomen_attitude(t,WK_R.f,WK_R); %time-varying 
% [R W W_dot theta_B] = body_attitude(15.65*pi/180); % fixed 
% [Q_A W_A W_A_dot theta_A] = abdomen_attitude(17.32*pi/180); % fixed

[L_R,L_L,D_R,D_L,M_R,M_L,F_rot_R,F_rot_L,M_rot_R,M_rot_L,alpha_R,alpha_L,U_alpha_R_dot,U_alpha_L_dot,U_R,U_L] = ...
    wing_QS_aerodynamics_alab(SIM,INSECT,W_R,W_L,W_R_dot,W_L_dot,...
    x_dot,R,W,Q_R,Q_L,Euler_R,Euler_R_dot,Euler_L,Euler_L_dot);

if SIM.cfdOnFlag == 1

    fprintf('### Inside eom ###\n')
    disp(['Changing to ',SIM.runPath])
    cd(SIM.runPath);
    if SIM.collectData == 0
        if t == 0
            disp('t == 0, timeStep = 1, inside Euler Method')
            system('rm fort2matlab_it.*dat');
            system('rm matlab2fort_it.*dat');
            fprintf('---> Matlab time = 0: Starting Stream\n')
            fprintf('---> t = %.12f \n',t)
            fprintf('---> timeStep = %i \n',timeStep)
            pause
            system('./cleanRun.sh');
            system(['./stream ',ctrlName,' > run.log &']);
            pause(10)
            dlmwrite('matlab2fort_it.dat',[x(1), x(2), x(3), reshape(Euler_R,1,3), WK_R.beta, theta_B]);
        else
            disp('---> Beyond first time step')
            dlmwrite('matlab2fort_it.dat',[x(1), x(2), x(3), reshape(Euler_R,1,3), WK_R.beta, theta_B]);
        end
        count=0;
        pause(1);
        while exist('fort2matlab_it.dat','file')==0
            disp('---> Waiting on CFD simulation')
            fprintf('---> t = %.12f \n',t)
            fprintf('---> timeStep = %i \n',timeStep)
            pause(2)
            count = count+1;
            if count == 10000000
                disp('no file was generated in allowed time');
                break
            end
        end
    end
    
    [F_I_CFD,F_B_CFD,F_R_CFD,M_I_CFD,M_B_CFD,M_R_CFD] = integrate_surf_CFD('./',timeStep,Q_R,R,x);
    
    F_I_CFD(1) = F_I_CFD(1) + F_I_CFD(1);
    F_I_CFD(2) = F_I_CFD(2) - F_I_CFD(2);
    F_I_CFD(3) = F_I_CFD(3) + F_I_CFD(3);
    
    F_B_CFD(1) = F_B_CFD(1) + F_B_CFD(1);
    F_B_CFD(2) = F_B_CFD(2) - F_B_CFD(2);
    F_B_CFD(3) = F_B_CFD(3) + F_B_CFD(3);
    
    system('rm matlab2fort_it.dat');
    system('rm fort2matlab_it.dat');
    cd ../../matlab
else
    F_I_CFD = [0 0 0]';
    F_B_CFD = [0 0 0]';
    F_R_CFD = [0 0 0]';
    M_I_CFD = [0 0 0]';
    M_B_CFD = [0 0 0]';
    M_R_CFD = [0 0 0]';
    M_L_CFD = [0 0 0]';
    
end

% Forces in respective (left and right) wing frames:
F_R=L_R+D_R+F_rot_R;
F_L=L_L+D_L+F_rot_L;
M_R=M_R+M_rot_R;
M_L=M_L+M_rot_L;

% F_R=L_R+D_R;
% F_L=L_L+D_L;
% M_R=M_R;
% M_L=M_L;
M_A=zeros(3,1);



F_R_CFD=F_R_CFD;
F_L_CFD=F_R_CFD; F_L_CFD(2) = -F_L_CFD(2);
M_R_CFD=M_R_CFD;
M_L_CFD=-M_R_CFD; M_L_CFD(2)=-M_L_CFD(2);
M_A_CFD=zeros(3,1);

% Set forces to QS forces
f_a=[R*Q_R*F_R + R*Q_L*F_L;
    hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R;
    M_L;
    M_A];

% Set forces to CFD forces
% f_a=[R*Q_R*F_R_CFD + R*Q_L*F_L_CFD;
%     hat(INSECT.mu_R)*Q_R*F_R_CFD + hat(INSECT.mu_L)*Q_L*F_L_CFD;
%     M_R_CFD;
%     M_L_CFD;
%     M_A_CFD];

f_a_1=f_a(1:3);
f_a_2=f_a(4:15);

% gravitational force and moment
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:3);
f_g_2=f_g(4:15);

% Euler-Lagrange equation
xi_1=x_dot;
xi_2=[W; W_R; W_L; W_A];
xi_2_dot=[W_dot; W_R_dot; W_L_dot; W_A_dot];

[JJ,KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

[JJ_11,JJ_12,JJ_21,JJ_22] = inertia_sub_decompose_3_12(JJ);
[LL_11,LL_12,LL_21,LL_22] = inertia_sub_decompose_3_12(LL);
[~, ~, ~, co_ad_22] = inertia_sub_decompose_3_12(co_ad);

xi_1_dot = JJ_11\( -JJ_12*xi_2_dot -LL_11*xi_1 - LL_12*xi_2 + f_a_1 + f_g_1);

f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1 + JJ_22*xi_2) ...
    + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
f_tau = [zeros(3,1); f_tau_2];
tau = blkdiag(zeros(3), Q_R, Q_L, Q_A)*f_tau_2;

% xi=[xi_1;xi_2];
% xi_dot=JJ\( co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau);
% disp(norm(xi_dot - [xi_1_dot; xi_2_dot]));

X_dot=[xi_1; xi_1_dot];

end