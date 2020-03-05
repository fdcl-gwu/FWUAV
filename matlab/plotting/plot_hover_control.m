%% Figures for controlled hover
addpath('../modules', '../sim_data', '../');
% load('sim_QS_x_hover_control.mat');

h_x=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
    hold on;
    plot(t*WK.f,des.x_fit_t(ii, :), 'k');
    patch_downstroke(h_x,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$x$','interpreter','latex');
% print(h_x, 'hover_control_pos', '-depsc');

h_x_dot=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
    hold on;
    plot(t*WK.f,des.x_dot_fit_t(ii, :), 'k');
    patch_downstroke(h_x_dot,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\dot x$','interpreter','latex');
% print(h_x_dot, 'hover_control_vel', '-depsc');

h_err = figure;
subplot(3,1,1);
plot(t*WK.f, pos_err(1,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
subplot(3,1,2);
plot(t*WK.f, pos_err(2,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
ylabel('$\Delta x$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f, pos_err(3,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
