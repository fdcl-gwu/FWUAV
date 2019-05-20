clear all;
close all;

WK.f=10.1;
WK.t_shift = 0.0408;
%WK.beta=0*pi/180;
N=5001;
T=5/WK.f;
t=linspace(0,T,N);

% WK.phi_m=50*pi/180;
% WK.phi_K=0.4;
% WK.phi_0=10*pi/180;
% 
% WK.theta_m=45*pi/180;
% WK.theta_C=5;
% WK.theta_0=0;
% WK.theta_a=0.3;
% 
% WK.psi_m=10*pi/180;
% WK.psi_N=2;
% WK.psi_a=0;
% WK.psi_0=0;


WK.type='Monarch';

for k=1:N
    E=wing_kinematics(t(k),WK);
    phi(k)=E(1);
    theta(k)=E(2);
    psi(k)=E(3);
end
    
figure;
subplot(3,1,1);
plot(t/T,phi*180/pi);
grid on;
ylabel('$\phi$','interpreter','latex');
%text(0.2,-40,'upstroke','fontname','times');
%text(0.7,-40,'downstroke','fontname','times');

set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);

subplot(3,1,2);
plot(t/T,theta*180/pi);
grid on;
ylabel('$\theta$','interpreter','latex');
set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);

subplot(3,1,3);
plot(t/T,psi*180/pi);
grid on;
ylabel('$\psi$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);

print('wing_kinematics_Monarch','-depsc');
!mv wing_kinematics.eps ../doc/Figs