function sim_QS_all_flex
% simulate the complete dynamics (x,R,Q_R,Q_L,Q_A) for given torque acting
% on the joint along with passive wing pitch motion

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='sim_QS_all';

load('./sim_data/other_insects/morp_MONARCH.mat', 'INSECT');
INSECT.Ktorsion = 1e-4; % Torsional spring coefficient

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
    U_t_R = torsion_terms(INSECT, Q_R(:,:,k), 0);
    U_t_L = torsion_terms(INSECT, Q_L(:,:,k), 0);
    E(k)=1/2* xi'*JJ* xi + U + U_t_R + U_t_L;
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
[~, tau_t_R] = torsion_terms(INSECT, Q_R, 0);
[~, tau_t_L] = torsion_terms(INSECT, Q_L, 0);
dU(7:9) = dU(7:9) - tau_t_R;
dU(10:12) = dU(10:12) - tau_t_L;

co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

xi_dot=JJ\(-KK*xi + co_ad*JJ*xi + 0.5*KK'*xi - dU);

R_dot = R*hat(W);
Q_R_dot = Q_R*hat(W_R);
Q_L_dot = Q_L*hat(W_L);
Q_A_dot = Q_A*hat(W_A);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_R_dot,9,1); reshape(Q_L_dot,9,1); reshape(Q_A_dot,9,1); xi_dot];
end
