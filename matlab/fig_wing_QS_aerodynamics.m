function fig_wing_QS_aerodynamics
close all;

WK.f=10.1;
WK.beta=30*pi/180;

WK.phi_m=60*pi/180;
WK.phi_K=0.9;
WK.phi_0=10*pi/180;

WK.theta_m=40*pi/180;
WK.theta_C=4;
WK.theta_0=20*pi/180;
WK.theta_a=0;

WK.psi_m=30*pi/180;
WK.psi_N=1;
WK.psi_a=0;
WK.psi_0=10*pi/180;

WK.type='BermanWang';
myfig(WK,'QS_sym',true);
WK.theta_a=0.3;
myfig(WK,'QS_adv',true);
WK.theta_a=-0.3;
myfig(WK,'QS_dly',true);

WK.type='Monarch';
WK.t_shift = 0.0408;
myfig(WK,'QS_Monarch',true);

end

function myfig(WK,filename,bool_print)
N=1001;
T=1/WK.f;
t=linspace(0,T,N);

load('morp_MONARCH');

e2=[0 1 0]';
x_dot=[1.2 0 0]';
R=expm(10*pi/180*hat(e2));
W=[0 0 0]';

for k=1:N
    [E_R(:,k) E_R_dot(:,k), E_R_ddot(:,k)]=wing_kinematics(t(k),WK);
    [E_L(:,k) E_L_dot(:,k), E_L_ddot(:,k)]=wing_kinematics(t(k),WK);        
    
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k) W_R_dot(:,k) W_L_dot(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_L(:,k), E_R_dot(:,k), E_L_dot(:,k), E_R_ddot(:,k), E_L_ddot(:,k));
    [L_R(:,k) L_L(:,k) D_R(:,k) D_L(:,k) M_R(:,k) M_L(:,k) ...
        F_rot_R(:,k) F_rot_L(:,k) M_rot_R(:,k) M_rot_L(:,k) ...
        alpha_R(k) alpha_L(:,k) U_alpha_dot_R(:,k) U_alpha_dot_L(:,k) U_R(:,k) U_L(:,k)]...
        =wing_QS_aerodynamics(MONARCH, W_R(:,k), W_L(:,k), W_R_dot(:,k), W_L_dot(:,k), x_dot, R, W, Q_R(:,:,k), Q_L(:,:,k));
    norm_U(k)=norm(U_R(:,k));
    alpha_R_dot(k)=U_alpha_dot_R(k)/norm(U_R(:,k));
    alpha_L_dot(k)=U_alpha_dot_L(k)/norm(U_L(:,k));
    
    L_R_in_B(:,k) = Q_R(:,:,k)*L_R(:,k);
    D_R_in_B(:,k) = Q_R(:,:,k)*D_R(:,k);
    F_rot_R_in_B(:,k) = Q_R(:,:,k)*F_rot_R(:,k);
    L_L_in_B(:,k) = Q_L(:,:,k)*L_L(:,k);
    D_L_in_B(:,k) = Q_L(:,:,k)*D_L(:,k);
    F_rot_L_in_B(:,k) = Q_L(:,:,k)*F_rot_L(:,k);
end
L_in_B=L_R_in_B +L_L_in_B;
D_in_B=D_R_in_B +D_L_in_B;
F_rot_in_B=F_rot_R_in_B +F_rot_L_in_B;


hE=figure;
subplot(3,1,1);
plot(t/T,E_R(1,:)*180/pi);
grid on;set(gca,'XTick',[0 0.5 1]);
ylabel('$\phi$','interpreter','latex');
text(0.2,-40,'upstroke','fontname','times');
text(0.7,-40,'downstroke','fontname','times');

subplot(3,1,2);
plot(t/T,E_R(2,:)*180/pi);
grid on;set(gca,'XTick',[0 0.5 1]);
ylabel('$\theta$','interpreter','latex');

subplot(3,1,3);
plot(t/T,E_R(3,:)*180/pi);
grid on;set(gca,'XTick',[0 0.5 1]);
ylabel('$\psi$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');

h_alpha=figure;
plot(t/T,alpha_R,'b');
ylabel('$\alpha$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
grid on;set(gca,'XTick',[0 0.5 1]);

h_U_alpha_dot=figure;
plot(t/T,U_alpha_dot_R);
ylabel('$\dot \alpha \| U\|$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
grid on;set(gca,'XTick',[0 0.5 1]);


h_U=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t/T,U_R(ii,:),'r',t/T,U_L(ii,:),'b');
grid on;set(gca,'XTick',[0 0.5 1]);
    hold on;
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$U$','interpreter','latex');

h_F=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t/T,L_R(ii,:),'r',t/T,D_R(ii,:),'b',t/T,F_rot_R(ii,:),'m--');
    grid on;set(gca,'XTick',[0 0.5 1]);
    hold on;
    plot(t/T,L_R(ii,:)+D_R(ii,:)+F_rot_R(ii,:),'k','LineWidth',1.2);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
hl=legend({'$L_R$','$D_R$','$F_{\mathrm{rot}}$','$F_{\mathrm{total}}$'});
set(hl,'interpreter','latex');

h_FB=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t/T,L_in_B(ii,:),'r',t/T,D_in_B(ii,:),'b',t/T,F_rot_in_B(ii,:),'m--');
    hold on;
    plot(t/T,L_in_B(ii,:)+D_in_B(ii,:)+F_rot_in_B(ii,:),'k','LineWidth',1.2);
    grid on;set(gca,'XTick',[0 0.5 1]);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
hl=legend({'$L_R$','$D_R$','$F_{\mathrm{rot}}$','$F_{\mathrm{total}}$'});
set(hl,'interpreter','latex');

if bool_print
    print(hE,[filename '_E'],'-depsc');
    print(h_alpha,[filename '_alpha'],'-depsc');
    print(h_U_alpha_dot,[filename '_U_alpha_dot'],'-depsc');
    print(h_U, [filename '_U'],'-depsc');
    print(h_F, [filename '_F_W'],'-depsc');
    print(h_FB, [filename '_F_B'],'-depsc');    
    
    evalin('base',['!mv ' filename '*.eps ../doc/Figs']);
end

end

