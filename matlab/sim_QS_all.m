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

N_period = 1;
N_single = 50;
N_iters = 10;
N_per_iter = N_single / N_iters;
integrator = 'crgr4'; % ode45, ode45ct, crgr2, crgr3, crgr4

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
    
    case 'ode45ct'
    % Continuous ode45
    tic;
    [~, X] = ode45(@(t,X) eom(INSECT,t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    time_ode45ct = toc;

    case 'crgr2'
    % 2nd order Crouch Grossman
    tic;
    for k=1:N-1
        X(k+1, :) = crouch_grossman_2nd(INSECT, X(k, :)', dt);
    end
    time_crgr2 = toc;

    case 'crgr3'
    % 3rd order Crouch Grossman
%     a = [0, 0, 0;
%         -1/24, 0, 0;
%         161/24, -6, 0;];
%     b = [1, -2/3, 2/3];
    a = [0, 0, 0;
        3/4, 0, 0;
        119/216, 17/108, 0;];
    b = [13/51, -2/3, 24/17];
    tic;
    for k=1:N-1
        X(k+1, :) = crouch_grossman_3rd(INSECT, X(k, :)', dt, a, b);
    end
    time_crgr3 = toc;
    
    case 'crgr4'
    % 4th order Crouch grossman
    a = [0, 0, 0, 0, 0;
        0.8177227988124852, 0, 0, 0, 0;
        0.3199876375476427, 0.0659864263556022, 0, 0, 0;
        0.9214417194464946, 0.4997857776773573, -1.0969984448371582, 0, 0;
        0.3552358559023322, 0.2390958372307326, 1.3918565724203246, -1.1092979392113465, 0];
    b = [0.1370831520630755, -0.0183698531564020, 0.7397813985370780, -0.1907142565505889, 0.3322195591068374];
    c = [0, 0.8177227988124852, 0.3859740639032449, 0.3242290522866937, 0.8768903263420429];
    tic;
    for k=1:N-1
        X(k+1, :) = crouch_grossman_4th(INSECT, X(k, :)', dt, a, b);
    end
    time_crgr4 = toc;
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
plot(t*f,(E-E(1))./E);
ylabel('E');

figure;
plot(t*f,err_R);
ylabel('deviation of R');

figure;
plot(t*f,err_QR);
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

function Xout = crouch_grossman_2nd(INSECT, X, dt)
%%
x1=X(1:3);
R1=reshape(X(4:12),3,3);
Q_R1=reshape(X(13:21),3,3);
Q_L1=reshape(X(22:30),3,3);
Q_A1=reshape(X(31:39),3,3);
x_dot1=X(40:42);
W1=X(43:45);
W_R1=X(46:48);
W_L1=X(49:51);
W_A1=X(52:54);

xi1=X(40:54);
K1_xi = deriv_xi(INSECT, x1, R1, Q_R1, Q_L1, Q_A1, x_dot1, W1, W_R1, W_L1, W_A1);

x2 = x1 + dt * x_dot1;
R2 = R1 * expmhat(dt*W1); Q_R2 = Q_R1 * expmhat(dt*W_R1); Q_L2 = Q_L1 * expmhat(dt*W_L1); Q_A2 = Q_A1 * expmhat(dt*W_A1);
xi2 = xi1 + dt * K1_xi;
x_dot2 = xi2(1:3); W2 = xi2(4:6); W_R2 = xi2(7:9); W_L2 = xi2(10:12); W_A2 = xi2(13:15);
K2_xi = deriv_xi(INSECT, x2, R2, Q_R2, Q_L2, Q_A2, x_dot2, W2, W_R2, W_L2, W_A2);

xout = x1 + dt/2 * (x_dot1+x_dot2);
Rout = R1 * expmhat(dt/2*W1) * expmhat(dt/2*W2);
Q_Rout = Q_R1 * expmhat(dt/2*W_R1) * expmhat(dt/2*W_R2);
Q_Lout = Q_L1 * expmhat(dt/2*W_L1) * expmhat(dt/2*W_L2);
Q_Aout = Q_A1 * expmhat(dt/2*W_A1) * expmhat(dt/2*W_A2);
xiout = xi1 + dt/2 * (K1_xi + K2_xi);

Xout = [xout; reshape(Rout, 9, 1); reshape(Q_Rout, 9, 1); reshape(Q_Lout, 9, 1); reshape(Q_Aout, 9, 1); xiout];
end

function Xout = crouch_grossman_3rd(INSECT, X, dt, a, b)
%%
x1=X(1:3);
R1=reshape(X(4:12),3,3);
Q_R1=reshape(X(13:21),3,3);
Q_L1=reshape(X(22:30),3,3);
Q_A1=reshape(X(31:39),3,3);
x_dot1=X(40:42);
W1=X(43:45);
W_R1=X(46:48);
W_L1=X(49:51);
W_A1=X(52:54);

xi1=X(40:54);
K1_xi = deriv_xi(INSECT, x1, R1, Q_R1, Q_L1, Q_A1, x_dot1, W1, W_R1, W_L1, W_A1);

x2 = x1 + dt*a(2,1)* x_dot1;
R2 = R1 * expmhat(dt*a(2,1)*W1); Q_R2 = Q_R1 * expmhat(dt*a(2,1)*W_R1);
Q_L2 = Q_L1 * expmhat(dt*a(2,1)*W_L1); Q_A2 = Q_A1 * expmhat(dt*a(2,1)*W_A1);
xi2 = xi1 + dt*a(2,1)* K1_xi;
x_dot2 = xi2(1:3); W2 = xi2(4:6); W_R2 = xi2(7:9); W_L2 = xi2(10:12); W_A2 = xi2(13:15);
K2_xi = deriv_xi(INSECT, x2, R2, Q_R2, Q_L2, Q_A2, x_dot2, W2, W_R2, W_L2, W_A2);

x3 = x1 + dt* (a(3,1)* x_dot1 + a(3,2)* x_dot2);
R3 = R1 * expmhat(dt*a(3,1)*W1) * expmhat(dt*a(3,2)*W2);
Q_R3 = Q_R1 * expmhat(dt*a(3,1)*W_R1) * expmhat(dt*a(3,2)*W_R2);
Q_L3 = Q_L1 * expmhat(dt*a(3,1)*W_L1) * expmhat(dt*a(3,2)*W_L2);
Q_A3 = Q_A1 * expmhat(dt*a(3,1)*W_A1) * expmhat(dt*a(3,2)*W_A2);
xi3 = xi1 + dt* (a(3,1)* K1_xi + a(3,2)* K2_xi);
x_dot3 = xi3(1:3); W3 = xi3(4:6); W_R3 = xi3(7:9); W_L3 = xi3(10:12); W_A3 = xi3(13:15);
K3_xi = deriv_xi(INSECT, x3, R3, Q_R3, Q_L3, Q_A3, x_dot3, W3, W_R3, W_L3, W_A3);

xout = x1 + dt* (b(1)* x_dot1 + b(2)* x_dot2 + b(3)* x_dot3);
Rout = R1 * expmhat(dt*b(1)*W1) * expmhat(dt*b(2)*W2) * expmhat(dt*b(3)*W3);
Q_Rout = Q_R1 * expmhat(dt*b(1)*W_R1) * expmhat(dt*b(2)*W_R2) * expmhat(dt*b(3)*W_R3);
Q_Lout = Q_L1 * expmhat(dt*b(1)*W_L1) * expmhat(dt*b(2)*W_L2) * expmhat(dt*b(3)*W_L3);
Q_Aout = Q_A1 * expmhat(dt*b(1)*W_A1) * expmhat(dt*b(2)*W_A2) * expmhat(dt*b(3)*W_A3);
xiout = xi1 + dt* (b(1)* K1_xi + b(2)* K2_xi + b(3)* K3_xi);

Xout = [xout; reshape(Rout, 9, 1); reshape(Q_Rout, 9, 1); reshape(Q_Lout, 9, 1); reshape(Q_Aout, 9, 1); xiout];
end

function Xout = crouch_grossman_4th(INSECT, X, dt, a, b)
%%
x1=X(1:3);
R1=reshape(X(4:12),3,3);
Q_R1=reshape(X(13:21),3,3);
Q_L1=reshape(X(22:30),3,3);
Q_A1=reshape(X(31:39),3,3);
x_dot1=X(40:42);
W1=X(43:45);
W_R1=X(46:48);
W_L1=X(49:51);
W_A1=X(52:54);

xi1=X(40:54);
K1_xi = deriv_xi(INSECT, x1, R1, Q_R1, Q_L1, Q_A1, x_dot1, W1, W_R1, W_L1, W_A1);

x2 = x1 + dt*a(2,1)* x_dot1;
R2 = R1 * expmhat(dt*a(2,1)*W1); Q_R2 = Q_R1 * expmhat(dt*a(2,1)*W_R1);
Q_L2 = Q_L1 * expmhat(dt*a(2,1)*W_L1); Q_A2 = Q_A1 * expmhat(dt*a(2,1)*W_A1);
xi2 = xi1 + dt*a(2,1)* K1_xi;
x_dot2 = xi2(1:3); W2 = xi2(4:6); W_R2 = xi2(7:9); W_L2 = xi2(10:12); W_A2 = xi2(13:15);
K2_xi = deriv_xi(INSECT, x2, R2, Q_R2, Q_L2, Q_A2, x_dot2, W2, W_R2, W_L2, W_A2);

x3 = x1 + dt* (a(3,1)* x_dot1 + a(3,2)* x_dot2);
R3 = R1 * expmhat(dt*a(3,1)*W1) * expmhat(dt*a(3,2)*W2);
Q_R3 = Q_R1 * expmhat(dt*a(3,1)*W_R1) * expmhat(dt*a(3,2)*W_R2);
Q_L3 = Q_L1 * expmhat(dt*a(3,1)*W_L1) * expmhat(dt*a(3,2)*W_L2);
Q_A3 = Q_A1 * expmhat(dt*a(3,1)*W_A1) * expmhat(dt*a(3,2)*W_A2);
xi3 = xi1 + dt* (a(3,1)* K1_xi + a(3,2)* K2_xi);
x_dot3 = xi3(1:3); W3 = xi3(4:6); W_R3 = xi3(7:9); W_L3 = xi3(10:12); W_A3 = xi3(13:15);
K3_xi = deriv_xi(INSECT, x3, R3, Q_R3, Q_L3, Q_A3, x_dot3, W3, W_R3, W_L3, W_A3);

x4 = x1 + dt* (a(4,1)* x_dot1 + a(4,2)* x_dot2 + a(4,3)* x_dot3);
R4 = R1 * expmhat(dt*a(4,1)*W1) * expmhat(dt*a(4,2)*W2) * expmhat(dt*a(4,3)*W3);
Q_R4 = Q_R1 * expmhat(dt*a(4,1)*W_R1) * expmhat(dt*a(4,2)*W_R2) * expmhat(dt*a(4,3)*W_R3);
Q_L4 = Q_L1 * expmhat(dt*a(4,1)*W_L1) * expmhat(dt*a(4,2)*W_L2) * expmhat(dt*a(4,3)*W_L3);
Q_A4 = Q_A1 * expmhat(dt*a(4,1)*W_A1) * expmhat(dt*a(4,2)*W_A2) * expmhat(dt*a(4,3)*W_A3);
xi4 = xi1 + dt* (a(4,1)* K1_xi + a(4,2)* K2_xi + a(4,3)* K3_xi);
x_dot4 = xi4(1:3); W4 = xi4(4:6); W_R4 = xi4(7:9); W_L4 = xi4(10:12); W_A4 = xi4(13:15);
K4_xi = deriv_xi(INSECT, x4, R4, Q_R4, Q_L4, Q_A4, x_dot4, W4, W_R4, W_L4, W_A4);

x5 = x1 + dt* (a(5,1)* x_dot1 + a(5,2)* x_dot2 + a(5,3)* x_dot3 + a(5,4)* x_dot4);
R5 = R1 * expmhat(dt*a(5,1)*W1) * expmhat(dt*a(5,2)*W2) * expmhat(dt*a(5,3)*W3) * expmhat(dt*a(5,4)*W4);
Q_R5 = Q_R1 * expmhat(dt*a(5,1)*W_R1) * expmhat(dt*a(5,2)*W_R2) * expmhat(dt*a(5,3)*W_R3) * expmhat(dt*a(5,4)*W_R4);
Q_L5 = Q_L1 * expmhat(dt*a(5,1)*W_L1) * expmhat(dt*a(5,2)*W_L2) * expmhat(dt*a(5,3)*W_L3) * expmhat(dt*a(5,4)*W_L4);
Q_A5 = Q_A1 * expmhat(dt*a(5,1)*W_A1) * expmhat(dt*a(5,2)*W_A2) * expmhat(dt*a(5,3)*W_A3) * expmhat(dt*a(5,4)*W_A4);
xi5 = xi1 + dt* (a(5,1)* K1_xi + a(5,2)* K2_xi + a(5,3)* K3_xi + a(5,4)* K4_xi);
x_dot5 = xi5(1:3); W5 = xi5(4:6); W_R5 = xi5(7:9); W_L5 = xi5(10:12); W_A5 = xi5(13:15);
K5_xi = deriv_xi(INSECT, x5, R5, Q_R5, Q_L5, Q_A5, x_dot5, W5, W_R5, W_L5, W_A5);

xout = x1 + dt* (b(1)* x_dot1 + b(2)* x_dot2 + b(3)* x_dot3 + b(4)* x_dot4 + b(5)* x_dot5);
Rout = R1 * expmhat(dt*b(1)*W1) * expmhat(dt*b(2)*W2) * expmhat(dt*b(3)*W3) * expmhat(dt*b(4)*W4) * expmhat(dt*b(5)*W5);
Q_Rout = Q_R1 * expmhat(dt*b(1)*W_R1) * expmhat(dt*b(2)*W_R2) * expmhat(dt*b(3)*W_R3) * expmhat(dt*b(4)*W_R4) * expmhat(dt*b(5)*W_R5);
Q_Lout = Q_L1 * expmhat(dt*b(1)*W_L1) * expmhat(dt*b(2)*W_L2) * expmhat(dt*b(3)*W_L3) * expmhat(dt*b(4)*W_L4) * expmhat(dt*b(5)*W_L5);
Q_Aout = Q_A1 * expmhat(dt*b(1)*W_A1) * expmhat(dt*b(2)*W_A2) * expmhat(dt*b(3)*W_A3) * expmhat(dt*b(4)*W_A4) * expmhat(dt*b(5)*W_A5);
xiout = xi1 + dt* (b(1)* K1_xi + b(2)* K2_xi + b(3)* K3_xi + b(4)* K4_xi + b(5)* K5_xi);

Xout = [xout; reshape(Rout, 9, 1); reshape(Q_Rout, 9, 1); reshape(Q_Lout, 9, 1); reshape(Q_Aout, 9, 1); xiout];
end

function xi_dot = deriv_xi(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
%%
    xi=[x_dot; W; W_R; W_L; W_A];
    [JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
    [~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);

    co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

    xi_dot=JJ\(-KK*xi + co_ad*JJ*xi + 0.5*KK'*xi - dU);
end
