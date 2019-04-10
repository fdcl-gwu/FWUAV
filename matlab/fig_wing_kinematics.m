clear all;
close all;

WK.f=50;
WK.beta=5*pi/180;
N=301;
T=1/WK.f;
t=linspace(0,T,N);

WK.phi_m=50*pi/180;
WK.phi_K=0.4;
WK.phi_0=10*pi/180;

WK.theta_m=45*pi/180;
WK.theta_C=5;
WK.theta_0=0;
WK.theta_a=0.3;

WK.psi_m=10*pi/180;
WK.psi_N=2;
WK.psi_a=0;
WK.psi_0=0;

i=0;
for phi_K=[0.001 1]
    WK.phi_K=phi_K;
    i=i+1;
    for k=1:N
        E=wing_kinematics(t(k),WK);
        phi(i,k)=E(1);
    end    
end

i=0;
for theta_C=[0.1 3 100]
    WK.theta_C=theta_C;
    i=i+1;
    for k=1:N
        E=wing_kinematics(t(k),WK);
        theta(i,k)=E(2);
    end    
end

i=0;
for psi_N=[1 2]
    WK.psi_N=psi_N;
    i=i+1;
    for k=1:N
        E=wing_kinematics(t(k),WK);
        psi(i,k)=E(3);
    end    
end
    
figure;
subplot(3,1,1);
plot(t/T,phi*180/pi);
grid on;
ylabel('$\phi$','interpreter','latex');
h_L=legend('$\phi_K=0.001$','$\phi_K=1$');
set(h_L,'interpreter','latex');
text(0.2,-40,'upstroke','fontname','times');
text(0.7,-40,'downstroke','fontname','times');

set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);

subplot(3,1,2);
plot(t/T,theta*180/pi);
grid on;
ylabel('$\theta$','interpreter','latex');
h_L=legend('$\theta_C=0.1$','$\theta_C=3$','$\theta_C=100$');
set(h_L,'interpreter','latex');
set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);
subplot(3,1,3);
plot(t/T,psi*180/pi);
grid on;
ylabel('$\theta$','interpreter','latex');
h_L=legend('$\psi_N=1$','$\psi_N=2$');
set(h_L,'interpreter','latex');
set(gca,'XTick',[0 0.5 1]);
set(gca,'FontSize',14);

print('wing_kinematics','-depsc');
!mv wing_kinematics.eps ../doc/Figs