function sim_QS_all
%% simulate the complete dynamics (x,R,Q_R,Q_L,Q_A).

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='sim_QS_all';

load('morp_MONARCH.mat', 'MONARCH');
INSECT=MONARCH;
f=10.2247;

rng default;
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

N_period = 50;
N_single = 100;
N_iters = 10;
integrator = 'ode45'; % ode45, crgr

N_per_iter = N_single / N_iters;
N = 1 + N_period * N_single;
t=linspace(0,1/f*N_period,N);
X = zeros(N, length(X0));
X(1, :) = X0;
dt = t(2) - t(1);

switch integrator
    case 'ode45'
    % ODE45
    tic;
    for con=1:(N_iters*N_period)
        idx = (1+(con-1)*N_per_iter):(1+con*N_per_iter);
        [~, X(idx, :)] = ode45(@(t,X) eom(INSECT,t,X), t(idx), X(idx(1),:)', ...
            odeset('AbsTol',1e-6,'RelTol',1e-6));
    end
%     [~, X] = ode45(@(t,X) eom(INSECT,t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    time_ode45 = toc;

    case 'crgr'
    % 2nd order Crouch Grossman
    tic;
    for k=1:N-1
        X(k+1, :) = crouch_grossman_X(INSECT, X(k, :)', dt);
    end
    time_crgr = toc;
end

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

I = eye(3);
for k=1:N
    xi=[x_dot(:,k); W(:,k); W_R(:,k); W_L(:,k); W_A(:,k)];    
    JJ = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k), W_A(:,k));
    U = potential(INSECT,x(:,k),R(:,:,k),Q_R(:,:,k),Q_L(:,:,k),Q_A(:,:,k));
    E(k)=1/2* xi'*JJ* xi + U;
    err_R(k) = norm(R(:,:,k)'*R(:,:,k) - I, 'fro');
    err_QR(k) = norm(Q_R(:,:,k)'*Q_R(:,:,k) - I, 'fro');
end

figure;
plot(t,(E-E(1))./E);
ylabel('E');

figure;
plot(t,err_R);
ylabel('deviation of R');

figure;
plot(t,err_QR);
ylabel('deviation of Q_R');

save(filename);
evalin('base',['load ' filename]);
end

function X_dot = eom(INSECT, t, X)
%%
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

xi=[x_dot; W; W_R; W_L; W_A];
[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);

co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

xi_dot=JJ\(-KK*xi + co_ad*JJ*xi + 0.5*KK'*xi - dU);

R_dot = R*hat(W);
Q_R_dot = Q_R*hat(W_R);
Q_L_dot = Q_L*hat(W_L);
Q_A_dot = Q_A*hat(W_A);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_R_dot,9,1); reshape(Q_L_dot,9,1); reshape(Q_A_dot,9,1); xi_dot];
end

function Xout = crouch_grossman_X(INSECT, X, dt)
%%
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

x1 = x + dt * x_dot;
R1 = R * expmhat(dt*W); Q_R1 = Q_R * expmhat(dt*W_R); Q_L1 = Q_L * expmhat(dt*W_L); Q_A1 = Q_A * expmhat(dt*W_A);
K1_xi = deriv_xi(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
xi1 = X(40:54) + dt * K1_xi;
x_dot1 = xi1(1:3); W1 = xi1(4:6); W_R1 = xi1(7:9); W_L1 = xi1(10:12); W_A1 = xi1(13:15);

x2 = x + dt/2 * (x_dot+x_dot1);
R2 = R * expmhat(dt/2*W) * expmhat(dt/2*W1);
Q_R2 = Q_R * expmhat(dt/2*W_R) * expmhat(dt/2*W_R1);
Q_L2 = Q_L * expmhat(dt/2*W_L) * expmhat(dt/2*W_L1);
Q_A2 = Q_A * expmhat(dt/2*W_A) * expmhat(dt/2*W_A1);
K2_xi = deriv_xi(INSECT, x1, R1, Q_R1, Q_L1, Q_A1, x_dot1, W1, W_R1, W_L1, W_A1);
xi2 = X(40:54) + dt/2 * (K1_xi + K2_xi);

Xout = [x2; reshape(R2, 9, 1); reshape(Q_R2, 9, 1); reshape(Q_L2, 9, 1); reshape(Q_A2, 9, 1); xi2];
end

function xi_dot = deriv_xi(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
%%
    xi=[x_dot; W; W_R; W_L; W_A];
    [JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
    [~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);

    co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

    xi_dot=JJ\(-KK*xi + co_ad*JJ*xi + 0.5*KK'*xi - dU);
end
