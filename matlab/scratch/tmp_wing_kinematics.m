%% verify time-derivatives of the functions "wing_kinematics" and "wing_attitude"
clear all;
close all;

%%
WK.type='BermanWang';
WK.f=50;
WK.beta=rand;
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

%%
WK.f=10.1;
WK.type='Monarch';
WK.t_shift=2.87128e-03;

N=10001;
T=5/WK.f;
t=linspace(0,T,N);

for k=1:N
    [E_R(:,k) E_R_dot(:,k) E_R_ddot(:,k)]=wing_kinematics(t(k),WK);
    [Q_R(:,:,k) Q_L(:,:,k) W_R(:,k) W_L(:,k) W_R_dot(:,k) W_L_dot(:,k)]=wing_attitude(WK.beta, E_R(:,k), E_R(:,k), E_R_dot(:,k), E_R_dot(:,k), E_R_ddot(:,k), E_R_ddot(:,k));
end

for k=1:N-1
    W_R_diff(:,k)=vee(Q_R(:,:,k)'*Q_R(:,:,k+1)-eye(3))/(t(k+1)-t(k));
    W_L_diff(:,k)=vee(Q_L(:,:,k)'*Q_L(:,:,k+1)-eye(3))/(t(k+1)-t(k));
end

figure;
my_ylabel={'$\dot\phi$','$\dot\theta$','$\dot\psi$'};
for ii=1:3
    subplot(3,1,ii);
    plot(t,E_R_dot(ii,:),'r',t(2:end),diff(E_R(ii,:))./diff(t),'b--');
    ylabel(my_ylabel{ii},'interpreter','latex');
end

figure;
my_ylabel={'$\ddot\phi$','$\ddot\theta$','$\ddot\psi$'};
for ii=1:3
    subplot(3,1,ii);
    plot(t,E_R_ddot(ii,:),'r',t(2:end),diff(E_R_dot(ii,:))./diff(t),'b--');
    ylabel(my_ylabel{ii},'interpreter','latex');
end

figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,W_R(ii,:),'r',t(2:end),W_R_diff(ii,:),'b--');
end
subplot(3,1,2);ylabel('$\Omega_R$','interpreter','latex');

figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,W_L(ii,:),'r',t(2:end),W_L_diff(ii,:),'b--');
end
subplot(3,1,2);ylabel('$\Omega_L$','interpreter','latex');

figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,W_R_dot(ii,:),'r',t(2:end),diff(W_R(ii,:)')'./diff(t),'b--');
end
subplot(3,1,2);ylabel('$\dot\Omega_R$','interpreter','latex');

figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,W_L_dot(ii,:),'r',t(2:end),diff(W_L(ii,:)')'./diff(t),'b--');
end
subplot(3,1,2);ylabel('$\dot\Omega_L$','interpreter','latex');
