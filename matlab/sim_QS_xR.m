function sim_QS_xR
% simulate the position and the attitude of thorax (x,R), 
% for given wing kinematics and abdomen attitude
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='sim_QS_xR';

load('morp_MONARCH');
INSECT=MONARCH;

WK.f=10.2;
WK.beta=30*pi/180;
WK.type='Monarch';
WK.ab_type='experimental';

N=501;
T=5/WK.f;
t=linspace(0,T,N);

x0=[0 0 0]';
R0=eye(3);
x_dot0=[1.0 0 -0.5]';

func_HH_rot_total = @(W02) max(abs(momentum(INSECT, WK, WK, 0, x0, R0, [x_dot0; [0; W02; 0]])));

tmpN=5001;
W02=linspace(-50,50,tmpN);
for k=1:tmpN
    HH_total_max(k)=func_HH_rot_total(W02(k));
end
plot(W02,HH_total_max);
[~, tmpI]=min(HH_total_max)
W02=W02(tmpI);

% options = optimoptions('fmincon');
% options.OptimalityTolerance = 1e-10;
% options.MaxFunctionEvaluations = 10000;
% options.MaxIterations = 10000; 
% fmincon(func_HH_rot_total, -20,[],[],[],[],[],[],[],options)


W0=[0 W02 0]';
X0=[x0; reshape(R0,9,1); x_dot0; W0];


[t X]=ode45(@(t,X) eom(INSECT, WK, WK, t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,13:15)';
W=X(:,16:18)';

R=zeros(3,3,N);
for k=1:N
    R(:,:,k)=reshape(X(k,4:12),3,3);
    [~, R(:,:,k) Q_R(:,:,k) Q_L(:,:,k) Q_A(:,:,k) theta_A(k) W(:,k) W_R(:,k) W_L(:,k) W_A(:,k) F_R(:,k) F_L(:,k) M_R(:,k) M_L(:,k) f_a(:,k) f_g(:,k) f_tau(:,k) tau(:,k)]= eom(INSECT, WK, WK, t(k), X(k,:)');
    F_B(:,k)=Q_R(:,:,k)*F_R(:,k) + Q_L(:,:,k)*F_L(:,k);
    
    [HH_rot_total(:,k) HH(:,k) ] = momentum(INSECT, WK, WK, t(k), x(:,k), R(:,:,k), [x_dot(:,k); W(:,k)]);
    [Euler_R(:,k), Euler_R_dot(:,k), Euler_R_ddot(:,k)] = wing_kinematics(t(k),WK);
end

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

end

function [HH_rot_total HH] = momentum(INSECT, WK_R, WK_L, t, x, R, xi_1)
x_dot=xi_1(1:3);
W=xi_1(4:6);

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R Q_L W_R W_L] = wing_attitude(WK_R.beta, Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A W_A] = abdomen_attitude(t,WK_R.f,WK_R);

xi_2=[W_R; W_L; W_A];
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);

HH = JJ*[xi_1; xi_2];
HH_rot_total = HH(4:6) + Q_R*HH(7:9) + Q_L*HH(10:12) + Q_A*HH(13:15);
end

function [X_dot R Q_R Q_L Q_A theta_A W W_R W_L W_A F_R F_L M_R M_L f_a f_g f_tau tau]= eom(INSECT, WK_R, WK_L, t, X)
x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);

% wing/abdoment attitude and aerodynamic force/moment
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R Q_L W_R W_L W_R_dot W_L_dot] = wing_attitude(WK_R.beta, Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A W_A W_A_dot theta_A] = abdomen_attitude(t,WK_R.f,WK_R);
%[Q_A W_A W_A_dot theta_A] = abdomen_attitude(30*pi/180);

[L_R L_L D_R D_L M_R M_L ...
    F_rot_R F_rot_L M_rot_R M_rot_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
F_R=L_R+D_R+F_rot_R;
F_L=L_L+D_L+F_rot_L;
M_R=M_R+M_rot_R;
M_L=M_L+M_rot_L;
M_R=zeros(3,1);
M_L=zeros(3,1);
F_A=zeros(3,1);
M_A=zeros(3,1);

f_a=[R*Q_R*F_R + R*Q_L*F_L;
    hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R;
    M_L;
    M_A];
f_a_1=f_a(1:6);
f_a_2=f_a(7:15);
f_a=zeros(15,1);

% gravitational force and moment
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:6);
f_g_2=f_g(7:15);
%f_g=zeros(15,1);

% Euler-Lagrange equation
xi_1=[x_dot; W]; 
xi_2=[W_R; W_L; W_A];
xi_2_dot=[W_R_dot; W_L_dot; W_A_dot];

[JJ KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

[JJ_11 JJ_12 JJ_21 JJ_22] = inertia_sub_decompose(JJ);
[LL_11 LL_12 LL_21 LL_22] = inertia_sub_decompose(LL);
[co_ad_11, ~, ~, co_ad_22] = inertia_sub_decompose(co_ad);

C=[zeros(3,9);
    -Q_R -Q_L -Q_A];

tmp_1 = -(co_ad_11*JJ_11-C*co_ad_22*JJ_21)*xi_1 + (LL_11-C*LL_21)*xi_1;
tmp_2 = -(JJ_12-C*JJ_22)*xi_2_dot + (co_ad_11*JJ_12-C*co_ad_22*JJ_22)*xi_2 ...
    -(LL_12-C*LL_22)*xi_2;
tmp_f = f_a_1+f_g_1 - C*(f_a_2+f_g_2);

xi_1_dot=(JJ_11-C*JJ_21)\(-tmp_1+tmp_2+tmp_f);

f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1+JJ_22*xi_2) ...
    + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
f_tau_1 = C*f_tau_2;
f_tau = [f_tau_1; f_tau_2];

tau = blkdiag(Q_R, Q_L, Q_A)*f_tau_2;

% xi=[xi_1;xi_2];
% xi_dot=JJ\( co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau);
% disp(norm(xi_dot - [xi_1_dot; xi_2_dot]));

R_dot = R*hat(W);
X_dot=[x_dot; reshape(R_dot,9,1); xi_1_dot];
end

function [JJ_11 JJ_12 JJ_21 JJ_22] = inertia_sub_decompose(JJ)
JJ_11 = JJ(1:6,1:6);
JJ_12 = JJ(1:6,7:15);
JJ_21 = JJ(7:15,1:6);
JJ_22 = JJ(7:15,7:15);
end
