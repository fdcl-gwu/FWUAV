function [JJ, KK] = inertia_wing_sub(m, mu, xi, J, R, Q, x_dot, W, W_i)
% Returns the inertia matrix and its derivative for the subsystem being
% considered

% R_dot=R*hat(W);
% Q_dot=Q*hat(W_i);
hat_mu = hat(mu);
hat_Q_xi = hat(Q*xi);
hat_xi = hat(xi);
hat_RT_x_dot = hat(R'*x_dot);
hat_W = hat(W);
hat_QT_W = hat(Q'*W);

JJ=zeros(9,9);

JJ(1:3,1:3)=m*eye(3);
JJ(1:3,4:6)=-m*R*(hat_mu+hat_Q_xi);
JJ(1:3,7:9)=-m*R*Q*hat_xi;

JJ(4:6,1:3)=JJ(1:3,4:6)';
JJ(4:6,4:6)=m*(hat_mu')*hat_mu+Q*J*Q'+m*(hat_mu'*hat_Q_xi+hat_Q_xi'*hat_mu);
JJ(4:6,7:9)=Q*J+m*hat_mu'*Q*hat_xi;

JJ(7:9,1:3)=JJ(1:3,7:9)';
JJ(7:9,4:6)=JJ(4:6,7:9)';
JJ(7:9,7:9)=J;

KK=zeros(9,9);

KK(1:3,4:6) = m*R*hat((hat_mu+hat_Q_xi)*W) + m*R*hat(Q*hat_xi*W_i);
KK(1:3,7:9) = -m*R*hat_W*Q*hat_xi + m*R*Q*hat(hat_xi*W_i);
KK(4:6,4:6) = m*(hat_mu+hat_Q_xi)*hat_RT_x_dot;
KK(4:6,7:9) = m*hat_RT_x_dot*Q*hat_xi - Q*hat(J*Q'*W) + Q*J*hat_QT_W ...
    -m*hat_mu*hat_W*Q*hat_xi - m* hat(hat_mu*W)*Q*hat_xi ...
    -Q*hat(J*W_i) + m*hat_mu*Q*hat(hat_xi*W_i);
KK(7:9,4:6) = m*hat_xi*Q'*hat_RT_x_dot;
KK(7:9,7:9) = m*hat_xi*hat(Q'*R'*x_dot) + J*hat_QT_W - m*hat_xi*hat(Q'*hat_mu*W);
end
