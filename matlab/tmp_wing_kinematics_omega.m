clear all;
close all;

syms phi theta psi real;
syms phi_dot theta_dot psi_dot;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

% %% right wing
% %phi=-rand;theta=-rand;psi=-rand;
% Q=[1 0 0;
%     0 cos(phi) -sin(phi);
%     0 sin(phi) cos(phi)]* ...
%     [cos(psi) sin(psi) 0;
%     -sin(psi) cos(psi) 0;
%     0 0 1]* ...
%     [cos(theta) 0 sin(theta);
%     0 1 0;
%     -sin(theta) 0 cos(theta)];
% 
% %Q_exp=expmso3(phi*e1)*expmso3(-psi*e3)*expmso3(theta*e2);
% %disp(norm(Q-Q_exp));
% 
% Q_dot=diff(Q,phi)*phi_dot + diff(Q,theta)*theta_dot + diff(Q,psi)*psi_dot;
% W=vee(simplify(Q.'*Q_dot));
% T=[diff(W,phi_dot), diff(W,theta_dot), diff(W,psi_dot)];
% simplify(W-T*[phi_dot; theta_dot; psi_dot])

%% left wing
%phi=-rand;theta=-rand;psi=-rand;
Q=[1 0 0;
    0 cos(phi) sin(phi);
    0 -sin(phi) cos(phi)]* ...
    [cos(psi) -sin(psi) 0;
    sin(psi) cos(psi) 0;
    0 0 1]* ...
    [cos(theta) 0 sin(theta);
    0 1 0;
    -sin(theta) 0 cos(theta)];

%Q_exp=expmso3(-phi*e1)*expmso3(psi*e3)*expmso3(theta*e2);
%disp(norm(Q-Q_exp));

Q_dot=diff(Q,phi)*phi_dot + diff(Q,theta)*theta_dot + diff(Q,psi)*psi_dot;
W=vee(simplify(Q.'*Q_dot));
T=[diff(W,phi_dot), diff(W,theta_dot), diff(W,psi_dot)];
simplify(W-T*[phi_dot; theta_dot; psi_dot])
