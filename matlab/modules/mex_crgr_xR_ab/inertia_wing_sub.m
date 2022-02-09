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
J_hat_QT_W = J*hat(Q'*W);
hat_mu_p_Q_xi = hat_mu+hat_Q_xi;
Q_H_Hxi_Wi = Q*hat(hat_xi*W_i);
Q_hat_xi = Q*hat_xi;

JJ=zeros(9,9);

JJ(1:3,1:3)=0.5*m*eye(3);
JJ(1:3,4:6)=-m*R*hat_mu_p_Q_xi;
JJ(1:3,7:9)=-m*R*Q_hat_xi;
JJ(4:6,4:6)=0.5*(m*(-hat_mu*hat_mu_p_Q_xi - hat_Q_xi*hat_mu) + Q*J*Q');
JJ(4:6,7:9)=Q*J-m*hat_mu*Q_hat_xi;
JJ(7:9,7:9)=0.5*J;

JJ = JJ+JJ';

KK=zeros(9,9);

KK(1:3,4:6) = m*R*hat(hat_mu_p_Q_xi*W + Q_hat_xi*W_i);
KK(1:3,7:9) = m*R*(-hat_W*Q_hat_xi + Q_H_Hxi_Wi);
KK(4:6,4:6) = m*hat_mu_p_Q_xi*hat_RT_x_dot;
KK(4:6,7:9) = m*(hat_RT_x_dot -hat_mu*hat_W - hat(hat_mu*W))*Q_hat_xi ...
    + Q*(-hat(J*(Q'*W + W_i)) + J_hat_QT_W) + m*hat_mu*Q_H_Hxi_Wi;
KK(7:9,4:6) = m*hat_xi*Q'*hat_RT_x_dot;
KK(7:9,7:9) = m*hat_xi*(hat(Q'*R'*x_dot) - hat(Q'*hat_mu*W)) + J_hat_QT_W;

end
