function fig_comp_ab
close all;

load sim_QS_x_wo_ab;
[pow, E, E_dot, eff] = compute_power(MONARCH.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);


h_x3=figure;
plot(x(1,:),x(3,:),'g');
hold on;

h_x=figure;
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'g');
    hold on;
end

h_v=figure;
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'g');
    hold on;
end

h_v3=figure;
plot(x_dot(1,:),x_dot(3,:),'g');
hold on;

h_theta=figure;
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'g');
hold on;
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'g');
hold on;

h_tau=figure;
subplot(3,1,1);
plot(t*WK.f, tau(4:6,:),'g');
hold on;
subplot(3,1,2);
plot(t*WK.f, tau(7:9,:),'g');
hold on;
subplot(3,1,3);
plot(t*WK.f, tau(11,:),'g');
hold on;

h_pow = figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,pow(ii,:),'g');
    hold on;
end

h_E = figure;
subplot(3,1,1);
plot(t*WK.f,E,'g');
hold on;
subplot(3,1,2);
plot(t*WK.f,E_dot,'g');
hold on;
subplot(3,1,3);
plot(t*WK.f,eff,'g');
hold on;


%%
load sim_QS_x_with_ab;
[pow, E, E_dot, eff] = compute_power(MONARCH.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);


figure(h_x3);
plot(x(1,:),x(3,:),'b');
set(gca,'YDir','reverse');
hold on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_3$','interpreter','latex');

figure(h_x);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'b');
end
xlabel('$t/T$','interpreter','latex');
subplot(2,1,1);
ylabel('$x_1$','interpreter','latex');
subplot(2,1,2);
ylabel('$x_3$','interpreter','latex');


figure(h_v);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'b');
end
xlabel('$t/T$','interpreter','latex');
subplot(2,1,1);
ylabel('$\dot x_1$','interpreter','latex');
subplot(2,1,2);
ylabel('$\dot x_3$','interpreter','latex');


figure(h_v3);
plot(x_dot(1,:),x_dot(3,:),'b');
set(gca,'YDir','reverse');
axis equal;

xlabel('$\dot x_1$','interpreter','latex');
ylabel('$\dot x_3$','interpreter','latex');

figure(h_theta);
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'b');
ylabel('$\theta_B$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'b');
ylabel('$\theta_A$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');

figure(h_tau);
subplot(3,1,1);
plot(t*WK.f, tau(4:6,:),'b');
ylabel('$\tau_R$','interpreter','latex');
subplot(3,1,2);
plot(t*WK.f, tau(7:9,:),'b');
ylabel('$\tau_L$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f, tau(11,:),'b');
ylabel('$\tau_A$','interpreter','latex');

figure(h_pow);
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,pow(ii,:),'b');
    hold on;
end

figure(h_E);
subplot(3,1,1);
plot(t*WK.f,E,'b');
subplot(3,1,2);
plot(t*WK.f,E_dot,'b');
subplot(3,1,3);
plot(t*WK.f,eff,'b');




%%
load sim_QS_x_with_opposite_ab;
[pow, E, E_dot, eff] = compute_power(MONARCH.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

figure(h_x3);
plot(x(1,:),x(3,:),'m');

figure(h_x);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'m');
end

figure(h_v);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'m');
    patch_downstroke(h_v,t*WK.f,Euler_R_dot);
    xlim([0 5]);
end

figure(h_v3);
plot(x_dot(1,:),x_dot(3,:),'m');
set(gca,'YDir','reverse');
axis equal;

figure(h_theta);
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'m');
xlim([0 5]);
patch_downstroke(h_theta,t*WK.f,Euler_R_dot);
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'m');
xlim([0 5]);
patch_downstroke(h_theta,t*WK.f,Euler_R_dot);

figure(h_tau);
subplot(3,1,1);
plot(t*WK.f, tau(4:6,:),'m');
xlim([0 5]);
patch_downstroke(h_tau,t*WK.f,Euler_R_dot);
subplot(3,1,2);
plot(t*WK.f, tau(7:9,:),'m');
xlim([0 5]);
patch_downstroke(h_tau,t*WK.f,Euler_R_dot);
subplot(3,1,3);
plot(t*WK.f, tau(11,:),'m');
xlim([0 5]);
patch_downstroke(h_tau,t*WK.f,Euler_R_dot);

figure(h_pow);
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,pow(ii,:),'m');
    xlim([0 5]);
    patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
    hold on;
end

figure(h_E);
subplot(3,1,1);
plot(t*WK.f,E,'m');
xlim([0 5]);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
subplot(3,1,2);
plot(t*WK.f,E_dot,'m');
xlim([0 5]);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
subplot(3,1,3);
plot(t*WK.f,eff,'m');
xlim([0 5]);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);


%%

bool_print=false;
figname='comp_ab';
if bool_print
    figure(h_x3);print([figname '_x3'],'-depsc');
    figure(h_x);print([figname '_x'],'-depsc');
    figure(h_v);print([figname '_v'],'-depsc');
    figure(h_v3);print([figname '_v3'],'-depsc');
    figure(h_theta);print([figname '_theta'],'-depsc');
    !mv comp_ab*.eps ../doc/Figs/
end

end


function [pow, E, E_dot, eff] = compute_power(m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A)

N=length(t);
pow=zeros(3,N);
for k=1:N
    pow(1,k) = tau(4:6,k)'*Q_R(:,:,k)*W_R(:,k); 
    pow(2,k) = tau(7:9,k)'*Q_L(:,:,k)*W_L(:,k); 
    pow(3,k) = tau(10:12,k)'*Q_A(:,:,k)*W_A(:,k); 
    E(k) = 0.5*m*x_dot(:,k)'*x_dot(:,k) - m*9.81*x(3);
end

E_dot = diff(E')./diff(t);
E_dot = [E_dot; E_dot(end)];

for k=1:N
    eff(k) = E_dot(k)/sum(pow(:,k));
end



end


