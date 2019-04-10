function [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK)
%wing_kinematics: compute wing Euler angles and their time-derivaties
%
% [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK) computes
%
%       Euler = [phi theta psi]' (flapping, pitch, deviation)
%
%for a given time and a struct variables with the following members
%describing wing kinematics
%
%         WK.f, Wk.beta
%         WK.phi_m, WK.phi_K, WK.phi_0
%         WK.theta_m, WK.theta_C, WK.theta_0, WK.theta_a
%         WK.psi_m, WK.psi_N, WK.psi_a, WK.psi_0
%

%% phi / flapping
A=WK.phi_m / asin(WK.phi_K);
a=WK.phi_K;
b=2*pi*WK.f;

phi = A*asin( a * cos(b*t)) + WK.phi_0;
phi_dot = -(A*a*b*sin(b*t))/(1 - a^2*cos(b*t)^2)^(1/2);
phi_ddot = (A*a*b^2*cos(b*t)*(a^2 - 1))/(1 - a^2*cos(b*t)^2)^(3/2);

%% theta / pitching
A=WK.theta_m / tanh(WK.theta_C);
a=WK.theta_C;
b=2*pi*WK.f;
c=WK.theta_a;

theta = A * tanh( a * sin(b*t + c) ) + WK.theta_0;
theta_dot = -A*a*b*cos(c + b*t)*(tanh(a*sin(c + b*t))^2 - 1);
theta_ddot = A*a*b^2*(tanh(a*sin(c + b*t))^2 - 1)*(sin(c + b*t) + 2*a*cos(c + b*t)^2*tanh(a*sin(c + b*t)));

%%  psi / deviation
A=WK.psi_m;
a=2*pi*WK.psi_N*WK.f;
b=WK.psi_a;

psi = A * cos( a*t + b ) + WK.psi_0;
psi_dot  = A * -a * sin(a*t+b);
psi_ddot = A * -a^2 * cos(a*t+b);

%% return values
Euler=[phi theta psi]';
Euler_dot=[phi_dot theta_dot psi_dot]';
Euler_ddot=[phi_ddot theta_ddot psi_ddot]';
end
