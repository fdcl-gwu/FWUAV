function tmp_inertia
close all;
filename='tmp_inertia';

m=rand;
mu=rand(3,1);
xi=rand(3,1);
J=rand(3,3);J=J+J';
R=expmso3(rand(3,1));
Q=expmso3(rand(3,1));
W=rand(3,1);
W_Q=rand(3,1);
x_dot=rand(3,1);

eps=1e-8;
[JJ JJ_dot KK]= inertia_wing_sub(m, mu, xi, J, R, Q, x_dot, W, W_Q);
[JJ_new ~]= inertia_wing_sub(m, mu, xi, J, R+eps*R*hat(W), Q+eps*Q*hat(W_Q), x_dot, W, W_Q);

JJ_dot-(JJ_new-JJ)/eps

JJ_dot*[x_dot; W; W_Q] - KK*[x_dot; W;W_Q]

save(filename);
evalin('base',['load ' filename]);
end

function [JJ_R JJ_L] = inertia_wing(INSECT, R, Q_R, Q_L)

JJ_R = inertia_wig_sub(INSECT.m_R, INSECT.mu_R, INSECT.xi_R, insect.J_R, R, Q_R);
JJ_L = inertia_wig_sub(INSECT.m_L, INSECT.mu_L, INSECT.xi_L, insect.J_L, R, Q_L);

end

function [JJ JJ_dot KK]= inertia_wing_sub(m, mu, xi, J, R, Q, x_dot, W, W_i)

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

JJ_dot=zeros(9,9);

JJ_dot(1:3,4:6)=-m*R_dot*(hat(mu)+hat(Q*xi))-m*R*(hat(Q_dot*xi));
JJ_dot(1:3,7:9)=-m*R_dot*Q*hat(xi)-m*R*Q_dot*hat(xi);

JJ_dot(4:6,1:3)=JJ_dot(1:3,4:6)';
JJ_dot(4:6,4:6)=Q_dot*J*Q' + Q*J*Q_dot' + m*(hat(mu)'*hat(Q_dot*xi)+hat(Q_dot*xi)'*hat(mu));
JJ_dot(4:6,7:9)=Q_dot*J + m*hat(mu)'*Q_dot*hat(xi);

JJ_dot(7:9,1:3)=JJ_dot(1:3,7:9)';
JJ_dot(7:9,4:6)=JJ_dot(4:6,7:9)';

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



