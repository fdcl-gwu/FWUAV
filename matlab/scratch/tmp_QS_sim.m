function tmp_QS_sim
global J J_i m_i m

J=rand*eye(3);
J_i=rand*eye(3);
m_i=rand;
m=rand;

R0=expmso3(rand(3,1));
Q_R0=expmso3(rand(3,1));
W0=rand(3,1);
W_R0=rand(3,1);
x_dot0=zeros(3,1);
x0=rand(3,1);

X0=[x0; reshape(R0,9,1); reshape(Q_R0,9,1); x_dot0; W0; W_R0];

f=10;
N=501;
t=linspace(0,1/f*5,N);
[t X]=ode45(@(t,X) eom(t,X), t, X0, odeset('AbsTol',1e-9,'RelTol',1e-9));

x=X(:,1:3)';
x_dot=X(:,22:24)';
W=X(:,25:27)';
W_R=X(:,28:30)';

R=zeros(3,3,N);Q_R=zeros(3,3,N);
for k=1:N
    R(:,:,k)=reshape(X(k,4:12),3,3);
    Q(:,:,k)=reshape(X(k,13:21),3,3);
end


for k=1:N
    xi=[x_dot(:,k); W(:,k); W_R(:,k)];
    JJ=zeros(3,3);
    JJ(1:3,1:3)=(m+m_i)*eye(3);
    JJ(4:9,4:9)=[Q(:,:,k)*J_i*Q(:,:,k)'+J Q(:,:,k)*J_i;
        J_i*Q(:,:,k)' J_i];       
    E(k)=1/2* xi'*JJ* xi;
end

figure;
plot(t,E);

end

function X_dot = eom(t, X)
global J J_i m m_i
x=X(1:3);
R=reshape(X(4:12),3,3);
Q=reshape(X(13:21),3,3);
x_dot=X(22:24);
W=X(25:27);
W_i=X(28:30);

JJ=zeros(3,3);
JJ(1:3,1:3)=(m+m_i)*eye(3);
JJ(4:9,4:9)=[Q*J_i*Q'+J Q*J_i;
    J_i*Q' J_i];

Kxi=[zeros(3,1);
    Q*hat(W_i)*J_i*Q'*W - Q*J_i*hat(W_i)*Q'*W + Q*hat(W_i)*J_i*W_i;
    -J_i*hat(W_i)*Q'*W];

dT=[zeros(3,1); zeros(3,1);
    -1/2*hat(J_i*Q'*W)*Q'*W - 1/2*hat(Q'*W)*J_i*Q'*W + hat(J_i*W_i)*Q'*W];

tmp=inv(JJ)*(-Kxi+dT);

R_dot=R*hat(W);
Q_dot=Q*hat(W_i);

X_dot=[x_dot; reshape(R_dot,9,1); reshape(Q_dot,9,1); tmp];
end












