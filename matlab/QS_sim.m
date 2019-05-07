function QS_sim
evalin('base','clear all');
close all;
filename='QS_sim';

INSECT.g=9.81;
INSECT.m_R=rand;
INSECT.mu_R=0.3*rand(3,1);
INSECT.xi_R=0.5*rand(3,1);
INSECT.J_R=rand_spd;
INSECT.m_L=rand;
INSECT.mu_L=0.3*rand(3,1);
INSECT.xi_L=0.5*rand(3,1);
INSECT.J_L=rand_spd;
INSECT.J=rand_spd;
INSECT.m=rand;

R0=expmso3(rand(3,1));
Q_R0=expmso3(rand(3,1));
Q_L0=expmso3(rand(3,1));
W0=rand(3,1);
W_R0=rand(3,1);
W_L0=rand(3,1);
x_dot0=zeros(3,1);
x0=rand(3,1);

X0=[x0; reshape(R0,9,1); reshape(Q_R0,9,1); reshape(Q_L0,9,1); x_dot0; W0; W_R0; W_L0];

f=10;
N=501;
t=linspace(0,1/f*5,N);
[t X]=ode45(@(t,X) eom(INSECT,t,X), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

x=X(:,1:3)';
x_dot=X(:,31:33)';
W=X(:,34:36)';
W_R=X(:,37:39)';
W_L=X(:,40:42)';

R=zeros(3,3,N);Q_R=zeros(3,3,N);
for k=1:N
    R(:,:,k)=reshape(X(k,4:12),3,3);
    Q_R(:,:,k)=reshape(X(k,13:21),3,3);
    Q_L(:,:,k)=reshape(X(k,22:30),3,3);
end

for k=1:N
    xi=[x_dot(:,k); W(:,k); W_R(:,k); W_L(:,k)];    
    JJ = inertia(INSECT, R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), x_dot(:,k), W(:,k), W_R(:,k), W_L(:,k));
    U = potential(INSECT,x(:,k),R(:,:,k),Q_R(:,:,k),Q_L(:,:,k));
    E(k)=1/2* xi'*JJ* xi + U;
end

figure;
plot(t,(E-E(1))./E);
    


save(filename);
evalin('base',['load ' filename]);
end

function X_dot = eom(INSECT, t, X)
x=X(1:3);
R=reshape(X(4:12),3,3);
Q_R=reshape(X(13:21),3,3);
Q_L=reshape(X(22:30),3,3);
x_dot=X(31:33);
W=X(34:36);
W_R=X(37:39);
W_L=X(40:42);

xi=[x_dot; W; W_R; W_L];
[JJ KK] = inertia(INSECT, R, Q_R, Q_L, x_dot, W, W_R, W_L);
[~, dU]=potential(INSECT,x,R,Q_R,Q_L);

co_ad=zeros(12,12);
co_ad(4:6,4:6) = -hat(W);
co_ad(7:9,7:9) = -hat(W_R);
co_ad(10:12,10:12) = -hat(W_L);

xi_dot=JJ\(-KK*xi + co_ad*xi + 0.5*KK'*xi - dU);

R_dot = R*hat(W);
Q_R_dot = Q_R*hat(W_R);
Q_L_dot = Q_L*hat(W_L);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_R_dot,9,1); reshape(Q_L_dot,9,1); xi_dot];
end

function [JJ KK] = inertia(INSECT, R, Q_R, Q_L, x_dot, W, W_R, W_L)
[JJ_R KK_R] = inertia_wing_sub(INSECT.m_R, INSECT.mu_R, INSECT.xi_R, INSECT.J_R, R, Q_R, x_dot, W, W_R);
[JJ_L KK_L] = inertia_wing_sub(INSECT.m_L, INSECT.mu_L, INSECT.xi_L, INSECT.J_L, R, Q_L, x_dot, W, W_L);

JJ=zeros(12,12);
JJ(1:3,1:3) = INSECT.m*eye(3) + JJ_R(1:3,1:3) + JJ_L(1:3,1:3);
JJ(1:3,4:6) = JJ_R(1:3,4:6) + JJ_L(1:3,4:6);
JJ(1:3,7:9) = JJ_R(1:3,7:9);
JJ(1:3,10:12) = JJ_L(1:3,7:9);
JJ(4:6,1:3) = JJ(1:3,4:6)';
JJ(4:6,4:6) = INSECT.J + JJ_R(4:6,4:6) + JJ_L(4:6,4:6);
JJ(4:6,7:9) = JJ_R(4:6,7:9);
JJ(4:6,10:12) = JJ_L(4:6,7:9);
JJ(7:9,1:3) = JJ(1:3,7:9)';
JJ(7:9,4:6) = JJ(4:6,7:9)';
JJ(7:9,7:9) = JJ_R(7:9,7:9);
JJ(7:9,10:12) = zeros(3);
JJ(10:12,1:3) = JJ(1:3,10:12)';
JJ(10:12,4:6) = JJ(4:6,10:12)';
JJ(10:12,7:9) = JJ(7:9,10:12)';
JJ(10:12,10:12) = JJ_L(7:9,7:9);

KK=zeros(12,12);
KK(1:3,4:6) = KK_R(1:3,4:6) + KK_L(1:3,4:6);
KK(1:3,7:9) = KK_R(1:3,7:9);
KK(1:3,10:12) = KK_L(1:3,7:9);
KK(4:6,4:6) = KK_R(4:6,4:6) + KK_L(4:6,4:6);
KK(4:6,7:9) = KK_R(4:6,7:9);
KK(4:6,10:12) = KK_L(4:6,7:9);
KK(7:9,4:6) = KK_R(7:9,4:6);
KK(7:9,7:9) = KK_R(7:9,7:9);
KK(10:12,4:6) = KK_L(7:9,4:6);
KK(10:12,10:12) = KK_L(7:9,7:9);
end

function [JJ KK] = inertia_wing_sub(m, mu, xi, J, R, Q, x_dot, W, W_i)

R_dot=R*hat(W);
Q_dot=Q*hat(W_i);

JJ=zeros(9,9);

JJ(1:3,1:3)=m*eye(3);
JJ(1:3,4:6)=-m*R*(hat(mu)+hat(Q*xi));
JJ(1:3,7:9)=-m*R*Q*hat(xi);

JJ(4:6,1:3)=JJ(1:3,4:6)';
JJ(4:6,4:6)=m*hat(mu)'*hat(mu)+Q*J*Q'+m*(hat(mu)'*hat(Q*xi)+hat(Q*xi)'*hat(mu));
JJ(4:6,7:9)=Q*J+m*hat(mu)'*Q*hat(xi);

JJ(7:9,1:3)=JJ(1:3,7:9)';
JJ(7:9,4:6)=JJ(4:6,7:9)';
JJ(7:9,7:9)=J;

KK=zeros(9,9);

KK(1:3,4:6) = m*R*hat((hat(mu)+hat(Q*xi))*W) + m*R*hat(Q*hat(xi)*W_i);
KK(1:3,7:9) = -m*R*hat(W)*Q*hat(xi) + m*R*Q*hat(hat(xi)*W_i);
KK(4:6,4:6) = m*(hat(mu)+hat(Q*xi))*hat(R'*x_dot);
KK(4:6,7:9) = m*hat(R'*x_dot)*Q*hat(xi) - Q*hat(J*Q'*W) + Q*J*hat(Q'*W) ...
    -m*hat(mu)*hat(W)*Q*hat(xi) - m* hat(hat(mu)*W)*Q*hat(xi) ...
    -Q*hat(J*W_i) + m*hat(mu)*Q*hat(hat(xi)*W_i);
KK(7:9,4:6) = m*hat(xi)*Q'*hat(R'*x_dot);
KK(7:9,7:9) = m*hat(xi)*hat(Q'*R'*x_dot) + J*hat(Q'*W) - m*hat(xi)*hat(Q'*hat(mu)*W);

end

function J=rand_spd

Q=expmso3(rand(3,1));
J=Q*diag(rand(3,1))*Q';
end

function [U dU]=potential(INSECT,x,R,Q_R,Q_L)
e3=[0 0 1]';

mg=INSECT.m*INSECT.g;
mg_R=INSECT.m_R*INSECT.g;
mg_L=INSECT.m_L*INSECT.g;
rho_R = INSECT.mu_R + Q_R*INSECT.xi_R;
rho_L = INSECT.mu_L + Q_L*INSECT.xi_L;

U_B = -mg*e3'*x;
U_R = -mg_R*e3' * (x + R*rho_R);
U_L = -mg_L*e3' * (x + R*rho_L);
U = U_B + U_R + U_L;

dU = [-(INSECT.m + INSECT.m_R + INSECT.m_L ) * INSECT.g * e3;
    mg_R*hat(R'*e3)*rho_R + mg_L*hat(R'*e3)*rho_L;
    mg_R*hat(Q_R'*R'*e3)*INSECT.xi_R;
    mg_L*hat(Q_L'*R'*e3)*INSECT.xi_L];

end





