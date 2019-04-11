%% verify time-derivatives of the functions "wing_kinematics" and "wing_attitude"
clear all;
close all;

WK.f=50;
WK.beta=rand;
N=10001;
T=1/WK.f;
t=linspace(0,T,N);

WK.phi_m=rand*pi/2;
WK.phi_K=rand;
WK.phi_0=rand*pi/2;

WK.theta_m=rand*pi/2;
WK.theta_C=rand;
WK.theta_0=rand;
WK.theta_a=rand;

WK.psi_m=0*pi/2;
WK.psi_N=2;
WK.psi_a=rand;
WK.psi_0=rand;

MONARCH.l=52.2E-3; %wing span
MONARCH.S=2.46E-3/2; %area of right wing
MONARCH.tilde_r_2 = 0.456*MONARCH.l; %non-dimensional radius of the second moment of wing area

for k=1:N
    [E_R(:,k) E_R_dot(:,k), E_R_ddot(:,k)]=wing_kinematics(t(k),WK);
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k) W_R_dot(:,k) W_L_dot(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_R(:,k), E_R_dot(:,k), E_R_dot(:,k), E_R_ddot(:,k), E_R_ddot(:,k));
    [alpha(k) alpha_dot(:,k) U_R(:,k) U_theta(k)]=wing_QS_aerodynamics(MONARCH, W_R(:,k), W_R_dot(:,k));
    F_rot(k)=alpha_dot(1,k)*norm(U_R(:,k));
    norm_U(k)=norm(U_R(:,k));
end

figure;
plot(t,alpha,'b',t,U_theta,'r--');
ylabel('$\alpha$','interpreter','latex');
figure;
plot(t,alpha_dot(1,:),'r', t,alpha_dot(2,:),'g' ,t(2:end),diff(alpha)./diff(t),'b--');
ylabel('$\dot\alpha$','interpreter','latex');

figure;
plot(U_R(1,:),U_R(3,:));
axis equal;
grid on;

figure;
plot(t,F_rot);
ylabel('$\dot \alpha \| U\|$','interpreter','latex');

figure;
plot(t,norm_U);
ylabel('$\| U\|$','interpreter','latex');


