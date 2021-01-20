clc; close all; clear all;


set(0, 'DefaultLineLineWidth', 2);

XLIM=[0 10];

load sim_QS_x;
t_QS_x=t; x_QS_x=X(1:3,:); xdot_QS_x=X(4:6,:); 
fb_QS_x=F_B_QS; fi_QS_x=F_I_QS; mr_QS_x=M_R_QS;
theta_B_QS_x = theta_B; Euler_R_dot_QS = Euler_R_dot;

load sim_QSR_x;
t_QSR_x=t; x_QSR_x=X(1:3,:); xdot_QSR_x=X(4:6,:);
fb_QSR_x=F_B_QS; fi_QSR_x=F_I_QS;  mr_QSR_x=M_R_QS;
theta_B_QSR_x = theta_B; Euler_R_dot_QSR = Euler_R_dot;

load 'sim_CFD_x';
M_L_CFD = -M_R; M_L(2,:)=-M_L(2,:);
t_CFD=t; x_CFD=X(1:3,:); xdot_CFD=X(4:6,:);
fb_CFD=F_B; fi_CFD=F_I; mr_CFD=M_R;

% Plot forces in body frame:
h_FB=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t_QS_x*WK.f,fb_QS_x(ii,:),'b-'); hold on;
    plot(t_QSR_x*WK.f,fb_QSR_x(ii,:),'r-'); hold on;
    plot(t_CFD*WK.f,fb_CFD(ii,:),'k-'); hold on;
    patch_downstroke(h_FB,t_QS_x*WK.f,Euler_R_dot_QS);
    xlim(XLIM);
end
h_FB;
subplot(3,1,1); ylabel('$F_{B,x}$ [N]','interpreter','latex','fontsize',10)
ylim([-0.02 0.02])
subplot(3,1,2); ylabel('$F_{B,y}$ [N]','interpreter','latex','fontsize',10)
legend('QS','QS w/ rot','NS','interpreter','latex','fontsize',10);
subplot(3,1,3); ylabel('$F_{B,z}$ [N]','interpreter','latex','fontsize',10)
subplot(3,1,3); xlabel('$t/T$','interpreter','latex','fontsize',10)
ylim([-0.03 0.025])

% Plot forces in inertial frame:
h_FI=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t_QS_x*WK.f,fi_QS_x(ii,:),'b-'); hold on;
    plot(t_QSR_x*WK.f,fi_QSR_x(ii,:),'r-'); hold on;
    plot(t_CFD*WK.f,fi_CFD(ii,:),'k-'); hold on;
    patch_downstroke(h_FI,t_QS_x*WK.f,Euler_R_dot_QS);
    xlim(XLIM);
end
h_FB;
subplot(3,1,1); ylabel('$F_{I,x}$ [N]','interpreter','latex','fontsize',10)
ylim([-0.02 0.025])
subplot(3,1,2); ylabel('$F_{I,y}$ [N]','interpreter','latex','fontsize',10)
legend('QS','QS w/ rot','NS','interpreter','latex','fontsize',10);
subplot(3,1,3); ylabel('$F_{I,z}$ [N]','interpreter','latex','fontsize',10)
subplot(3,1,3); xlabel('$t/T$','interpreter','latex','fontsize',10)
ylim([-0.03 0.025])

% Plot positions in inertial frame:
h_x=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t_QS_x*WK.f,x_QS_x(ii,:),'b-'); hold on
    plot(t_QSR_x*WK.f,x_QSR_x(ii,:),'r-'); hold on;
    plot(t_CFD*WK.f,x_CFD(ii,:),'k-'); hold on;
    patch_downstroke(h_x,t_QS_x*WK.f,Euler_R_dot_QS);
    xlim(XLIM);
end
subplot(3,1,1); ylabel('$x_{I}$ [m]','interpreter','latex','fontsize',10);
subplot(3,1,2); ylabel('$y_{I}$ [m]','interpreter','latex','fontsize',10)
legend('QS','QS w/ rot','NS','interpreter','latex','fontsize',10);
subplot(3,1,3); ylabel('$z_{I}$ [m]','interpreter','latex','fontsize',10)
subplot(3,1,3); xlabel('$t/T$','interpreter','latex','fontsize',10)

% Plot body pitch in inertial frame:
h_R=figure;
plot(t_QS_x*WK.f,theta_B_QS_x,'b-'); hold on
xlim(XLIM);
ylabel('$\theta_{B}$ [deg]','interpreter','latex','fontsize',10);
legend('QS','interpreter','latex','fontsize',10)
xlabel('$t/T$','interpreter','latex','fontsize',10)

h_MR=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t_QS_x*WK.f,mr_QS_x(ii,:),'b-'); hold on;
	plot(t_QS_x*WK.f,mr_QSR_x(ii,:),'r-'); hold on;
    plot(t_CFD*WK.f,mr_CFD(ii,:),'k-'); hold on;
    patch_downstroke(h_MR,t_QS_x*WK.f,Euler_R_dot_QS);
    xlim(XLIM);
end
h_MR;
subplot(3,1,1); ylabel('$M_{R,x}$ [N-m]','interpreter','latex','fontsize',10)
ylim([-1 1]*1e-3);
subplot(3,1,2); ylabel('$M_{R,y}$ [N-m]','interpreter','latex','fontsize',10)
ylim([-1 2]*1e-4);
legend('QS','QS w/ rot','NS','interpreter','latex','fontsize',10);
subplot(3,1,3); ylabel('$M_{R,z}$ [N-m]','interpreter','latex','fontsize',10)
subplot(3,1,3); xlabel('$t/T$','interpreter','latex','fontsize',10)
