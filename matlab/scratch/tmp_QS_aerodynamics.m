%% verify time-derivatives of the functions "wing_kinematics" and "wing_attitude"
addpath('../');

clear all;
close all;

% WK.f=50;
% WK.beta=rand;
% 
% WK.phi_m=rand*pi/2;
% WK.phi_K=rand;
% WK.phi_0=rand*pi/2;
% 
% WK.theta_m=rand*pi/2;
% WK.theta_C=rand;
% WK.theta_0=rand;
% WK.theta_a=rand;
% 
% WK.psi_m=rand*pi/2;
% WK.psi_N=2;
% WK.psi_a=rand;
% WK.psi_0=rand;

WK.f=50;
WK.beta=20*pi/180;

WK.phi_m=50*pi/180;
WK.phi_K=0.4;
WK.phi_0=10*pi/180;

WK.theta_m=45*pi/180;
WK.theta_C=2;
WK.theta_0=0;
WK.theta_a=0.3;

WK.psi_m=0*pi/180;
WK.psi_N=2;
WK.psi_a=0;
WK.psi_0=0;

N=10001;
T=1/WK.f;
t=linspace(0,T,N);

load('morp_MONARCH');

for k=1:N
    [E_R(:,k) E_R_dot(:,k), E_R_ddot(:,k)]=wing_kinematics(t(k),WK);
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k) W_R_dot(:,k) W_L_dot(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_R(:,k), E_R_dot(:,k), E_R_dot(:,k), E_R_ddot(:,k), E_R_ddot(:,k));
    [L_R(:,k) D_R(:,k) M_R(:,k) F_rot_R(:,k) M_rot_R(:,k) alpha(k) U_alpha_dot(:,k) U_R(:,k)]=wing_QS_aerodynamics(MONARCH, W_R(:,k), W_R_dot(:,k));
    norm_U(k)=norm(U_R(:,k));
    alpha_dot(k)=U_alpha_dot(k)/norm(U_R(:,k));
    
    L_R_in_B(:,k) = Q_R(:,:,k)*L_R(:,k);
    D_R_in_B(:,k) = Q_R(:,:,k)*D_R(:,k);
    F_rot_R_in_B(:,k) = Q_R(:,:,k)*F_rot_R(:,k);

end



figure;
plot(t,alpha,'b');
ylabel('$\alpha$','interpreter','latex');
figure;
plot(t,alpha_dot,'r');
hold on;
plot(t(2:end),diff(alpha)./diff(t),'b--');
ylabel('$\dot\alpha$','interpreter','latex');

figure;
plot(U_R(1,:),U_R(3,:));
axis equal;
grid on;

figure;
plot(t,U_alpha_dot);
ylabel('$\dot \alpha \| U\|$','interpreter','latex');

% figure;
% plot(t,norm_U);
% ylabel('$\| U\|$','interpreter','latex');

figure;

for ii=1:3
    subplot(3,1,ii);
    plot(t,L_R(ii,:),'r',t,D_R(ii,:),'b',t,F_rot_R(ii,:),'m--');
    hold on;
    plot(t,L_R(ii,:)+D_R(ii,:)+F_rot_R(ii,:),'k','LineWidth',1.2);
end
subplot(3,1,1);
hl=legend({'$L_R$','$D_R$','$F_{\mathrm{rot}}$','$F_{\mathrm{total}}$'});
set(hl,'interpreter','latex');
title('Aerodynamic forces in the wing-fixed frame');

figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,L_R_in_B(ii,:),'r',t,D_R_in_B(ii,:),'b',t,F_rot_R_in_B(ii,:),'m--');
    hold on;
    plot(t,L_R_in_B(ii,:)+D_R_in_B(ii,:)+F_rot_R_in_B(ii,:),'k','LineWidth',1.2);
end
subplot(3,1,1);
hl=legend({'$L_R$','$D_R$','$F_{\mathrm{rot}}$','$F_{\mathrm{total}}$'});
set(hl,'interpreter','latex');
title('Aerodynamic forces in the body-fixed frame');
