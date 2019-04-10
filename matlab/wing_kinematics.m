function [Euler, Euler_dot] = wing_kinematics(t,WK)
%wing_kinematics: compute wing Euler angles and their time-derivaties
%
% [Euler, Euler_dot] = wing_kinematics(t,WK) computes
%
%       Euler = [phi theta psi]' (flapping, pitch, deviation)
%       Euler_dot = [phi_dot theta_dot psi_dot]'
%
%for a given time and a struct variables with the following members
%describing wing kinematics
%
%         WK.f, Wk.beta
%         WK.phi_m, WK.phi_K, WK.phi_0
%         WK.theta_m, WK.theta_C, WK.theta_0, WK.theta_a
%         WK.psi_m, WK.psi_N, WK.psi_a, WK.psi_0
%

phi = WK.phi_m / asin(WK.phi_K) * asin( WK.phi_K * cos(2*pi*WK.f*t) ) + WK.phi_0;
theta = WK.theta_m / tanh(WK.theta_C) * tanh( WK.theta_C * sin(2*pi*WK.f*t + WK.theta_a) ) + WK.theta_0;
psi = WK.psi_m * cos( 2*pi*WK.psi_N*WK.f*t + WK.psi_a ) + WK.psi_0;

phi_dot = WK.phi_m / asin(WK.phi_K) * - (WK.phi_K * 2*pi*WK.f) * sin( 2*pi*WK.f*t) / sqrt(1 - WK.phi_K^2 * cos(2*pi*WK.f*t)^2);
theta_dot = WK.theta_m / tanh(WK.theta_C) * - WK.theta_C * 2*pi*WK.f *cos(2*pi*WK.f*t + WK.theta_a) * ...
    ( tanh( WK.theta_C * sin(2*pi*WK.f*t + WK.theta_a) )^2 -1);
psi_dot  = WK.psi_m * - 2*pi*WK.psi_N*WK.f * sin( 2*pi*WK.psi_N*WK.f*t + WK.psi_a );

Euler=[phi theta psi]';
Euler_dot=[phi_dot theta_dot psi_dot]';
end

