function sim_QS_all_flex
% simulate the complete dynamics (x,R,Q_R,Q_L,Q_A) for given torque acting
% on the joint along with passive wing pitch motion

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='sim_QS_all';

load('./sim_data/other_insects/morp_MONARCH.mat', 'INSECT');
INSECT.Kflex = 1e-4; % Torsional spring coefficient

R0=expmhat(rand(3,1));
Q_R0=expmhat(rand(3,1));
Q_L0=expmhat(rand(3,1));
Q_A0=expmhat(rand(3,1));
W0=rand(3,1);
W_R0=rand(3,1);
W_L0=rand(3,1);
W_A0=rand(3,1);
x_dot0=zeros(3,1);
x0=rand(3,1);

X0=[x0; reshape(R0,9,1); reshape(Q_R0,9,1); reshape(Q_L0,9,1); reshape(Q_A0,9,1);...
    x_dot0; W0; W_R0; W_L0; W_A0];

f=INSECT.f;
N=501;
N_period = 5;
t=linspace(0,1/f*N_period,N);
[t X]=ode45(@(t,X) eom(INSECT,t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,40:42)';
W=X(:,43:45)';
W_R=X(:,46:48)';
W_L=X(:,49:51)';
W_A=X(:,52:54)';

R=zeros(3,3,N);Q_R=zeros(3,3,N);
for k=1:N
    R(:,:,k)=reshape(X(k,4:12),3,3);
    Q_R(:,:,k)=reshape(X(k,13:21),3,3);
    Q_L(:,:,k)=reshape(X(k,22:30),3,3);
    Q_A(:,:,k)=reshape(X(k,31:39),3,3);
end

for k=1:N
    xi=[x_dot(:,k); W(:,k); W_R(:,k); W_L(:,k); W_A(:,k)];    
    JJ = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    U = potential(INSECT,x(:,k),R(:,:,k),Q_R(:,:,k),Q_L(:,:,k),Q_A(:,:,k));
    E(k)=1/2* xi'*JJ* xi + U;
end

figure;
plot(t*f,(E-E(1))./E);

save(filename);
evalin('base',['load ' filename]);
end

function X_dot = eom(INSECT, t, X)
x=X(1:3);
R=reshape(X(4:12),3,3);
Q_R=reshape(X(13:21),3,3);
Q_L=reshape(X(22:30),3,3);
Q_A=reshape(X(31:39),3,3);
x_dot=X(40:42);
W=X(43:45);
W_R=X(46:48);
W_L=X(49:51);
W_A=X(52:54);

xi = [x_dot; W; W_R; W_L; W_A];
[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
[~, dU] = potential(INSECT,x,R,Q_R,Q_L,Q_A);

co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

xi_dot=JJ\(-KK*xi + co_ad*JJ*xi + 0.5*KK'*xi - dU);

R_dot = R*hat(W);
Q_R_dot = Q_R*hat(W_R);
Q_L_dot = Q_L*hat(W_L);
Q_A_dot = Q_A*hat(W_A);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_R_dot,9,1); reshape(Q_L_dot,9,1); reshape(Q_A_dot,9,1); xi_dot];
end

function [U, dU] = potential(INSECT,x,R,Q_R,Q_L,Q_A)
% Returns the potential energy and the corresponding gravitational force
e3=[0 0 1]';

mg_B=INSECT.m_B*INSECT.g;
mg_R=INSECT.m_R*INSECT.g;
mg_L=INSECT.m_L*INSECT.g;
mg_A=INSECT.m_A*INSECT.g;

tmp_R = INSECT.mu_R + Q_R*INSECT.nu_R;
tmp_L = INSECT.mu_L + Q_L*INSECT.nu_L;
tmp_A = INSECT.mu_A + Q_A*INSECT.nu_A;

[U_flex_R, tau_flex_R] = flex_terms(INSECT, Q_R);
[U_flex_L, tau_flex_L] = flex_terms(INSECT, Q_L);

U_B = -mg_B*e3'*x;
U_R = -mg_R*e3' * (x + R*tmp_R);
U_L = -mg_L*e3' * (x + R*tmp_L);
U_A = -mg_A*e3' * (x + R*tmp_A);
U = U_B + U_R + U_L + U_A + U_flex_R + U_flex_L;

hat_RT_e3 = hat(R'*e3);
dU = [-(INSECT.m_B + INSECT.m_R + INSECT.m_L + INSECT.m_A) * INSECT.g * e3;
    mg_R*hat_RT_e3*tmp_R + mg_L*hat_RT_e3*tmp_L + mg_A*hat_RT_e3*tmp_A;
    mg_R*hat(Q_R'*R'*e3)*INSECT.nu_R - tau_flex_R;
    mg_L*hat(Q_L'*R'*e3)*INSECT.nu_L - tau_flex_L;
    mg_A*hat(Q_A'*R'*e3)*INSECT.nu_A];
end

function [U_flex, tau_flex, theta] = flex_terms(INSECT, Q)
beta = 0;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

sx = Q' * expmhat(beta*e2) * e1;
theta = atan2(e3'*sx, e1'*sx);
U_flex = 0.5 * INSECT.Kflex * theta^2;
% tau_flex = (INSECT.Kflex*theta)/((1+tan(theta)^2) * (e1'*sx)^2) * hat(sx)^2 * e2;
tau_flex = (INSECT.Kflex*theta)/((1+tan(theta)^2) * (e1'*sx)^2) * ((e2'*sx)*sx - e2);
end
