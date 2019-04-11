clear all;
close all;

syms A a b c t;
phi = A*asin(a*cos(b*t))
phi_dot = simplify(diff(phi,t))
phi_ddot = simplify(diff(phi_dot,t))

theta = A * tanh( a * sin(b*t + c) )
theta_dot = simplify(diff(theta,t))
theta_ddot = simplify(diff(theta_dot,t))
