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

for k=1:N
    [E_R(:,k) E_R_dot(:,k)]=wing_kinematics(t(k),WK);
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_R(:,k), E_R_dot(:,k), E_R_dot(:,k));
end

for k=1:N-1
    W_R_diff(:,k)=vee(Q_R(:,:,k)'*Q_R(:,:,k+1)-eye(3))/(t(k+1)-t(k));
    W_L_diff(:,k)=vee(Q_L(:,:,k)'*Q_L(:,:,k+1)-eye(3))/(t(k+1)-t(k));
end

figure;
subplot(3,1,1);
plot(t,E_R_dot(1,:),'r',t(2:end),diff(E_R(1,:))./diff(t),'b--');
ylabel('$\dot\phi$','interpreter','latex');
subplot(3,1,2);
plot(t,E_R_dot(2,:),'r',t(2:end),diff(E_R(2,:))./diff(t),'b--');
ylabel('$\dot\theta$','interpreter','latex');
subplot(3,1,3);
plot(t,E_R_dot(3,:),'r',t(2:end),diff(E_R(3,:))./diff(t),'b--');
ylabel('$\dot\psi$','interpreter','latex');

figure;
plot(t,W_R,'r',t(2:end),W_R_diff,'b--');
ylabel('$\Omega_R$','interpreter','latex');
figure;
plot(t,W_L,'r',t(2:end),W_L_diff,'b--');
ylabel('$\Omega_L$','interpreter','latex');