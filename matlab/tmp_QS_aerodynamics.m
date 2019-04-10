%% verify time-derivatives of the functions "wing_kinematics" and "wing_attitude"
clear all;
close all;

WK.f=50;
WK.beta=rand;
N=5001;
T=5/WK.f;
t=linspace(0,T,N);

WK.phi_m=rand;
WK.phi_K=rand;
WK.phi_0=rand;

WK.theta_m=rand;
WK.theta_C=rand;
WK.theta_0=rand;
WK.theta_a=rand;

WK.psi_m=rand;
WK.psi_N=2;
WK.psi_a=rand;
WK.psi_0=rand;

MONARCH.l=52.2E-3; %wing span
MONARCH.S=2.46E-3/2; %area of right wing
MONARCH.tilde_r_2 = 0.456*MONARCH.l; %non-dimensional radius of the second moment of wing area


for k=1:N
    [E_R(:,k) E_R_dot(:,k)]=wing_kinematics(t(k),WK);
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_R(:,k), E_R_dot(:,k), E_R_dot(:,k));
end
