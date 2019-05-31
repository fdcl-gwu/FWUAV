function fit_VICON_data
clear all;
close all;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

load VICON_data

N=length(t);
R = zeros(3,3,N);
R_R = zeros(3,3,N);
R_L = zeros(3,3,N);
R_A = zeros(3,3,N);
Q_R = zeros(3,3,N);
Q_L = zeros(3,3,N);
Q_A = zeros(3,3,N);

WR = 0.5*(T1+T2); % wing root

for k=1:N
    %% right wing
    % fit a plane passing thorough WR to RW1, RW2, RW3
    M_RR = (RW1(:,k)-WR(:,k))*(RW1(:,k)-WR(:,k))'...
        + (RW2(:,k)-WR(:,k))*(RW2(:,k)-WR(:,k))'...
        + (RW3(:,k)-WR(:,k))*(RW3(:,k)-WR(:,k))';
    [V D]=eig(M_RR);
    [~, I_min]=min(diag(D));
    n=V(:,I_min);
    tmp3 = cross(RW1(:,k)-WR(:,k), RW3(:,k)-WR(:,k));
    tmp3 = tmp3/norm(tmp3);
    if dot(tmp3, n) < 0
        n= -n;
    end
    
    R_R(:,3,k) = n;
    % project wing tip to the fitted plane
    tmp2 = (eye(3)-n*n')*(RW2(:,k)-WR(:,k));
    R_R(:,2,k) = tmp2/norm(tmp2);
    R_R(:,1,k) = cross(R_R(:,2,k),R_R(:,3,k));
    
    
    %% left wing
    % fit a plane passing thorough WR to LW1, LW2, LW3
    M_LL = (LW1(:,k)-WR(:,k))*(LW1(:,k)-WR(:,k))'...
        + (LW2(:,k)-WR(:,k))*(LW2(:,k)-WR(:,k))'...
        + (LW3(:,k)-WR(:,k))*(LW3(:,k)-WR(:,k))';
    [V D]=eig(M_LL);
    [~, I_min]=min(diag(D));
    n=V(:,I_min);
    tmp3 = cross(LW3(:,k)-WR(:,k), LW1(:,k)-WR(:,k));
    tmp3 = tmp3/norm(tmp3);
    if dot(tmp3, n) < 0
        n= -n;
    end
    
    R_L(:,3,k) = n;
    % project wing tip to the plane
    tmp2 = (eye(3)-n*n')*(WR(:,k)-LW2(:,k));
    R_L(:,2,k) = tmp2/norm(tmp2);
    R_L(:,1,k) = cross(R_L(:,2,k),R_L(:,3,k));
    
    %% thorax attitude
    R(:,1,k) = (T1(:,k) - T2(:,k))/norm(T1(:,k) - T2(:,k));
    R(:,2,k) = e2; % assume there is no body roll
    R(:,3,k) = cross( R(:,1,k), R(:,2,k));
    
    %% abdomen attitude
    tmpA = (A1(:,k) + A2(:,k))/2;
    R_A(:,1,k) = (T2(:,k) - tmpA)/norm(T2(:,k) - tmpA);
    R_A(:,2,k) = e2; % assume there is no body roll
    R_A(:,3,k) = cross( R_A(:,1,k), R_A(:,2,k));    
    
    %% relative attitude
    Q_R(:,:,k)= R(:,:,k)'*R_R(:,:,k);
    Q_L(:,:,k)= R(:,:,k)'*R_L(:,:,k);
    Q_A(:,:,k)= R(:,:,k)'*R_A(:,:,k);
    
    theta_th(k)=atan2(-R(3,1,k),R(1,1,k));
    theta_ab(k)=atan2(-Q_A(3,1,k),Q_A(1,1,k));
end

%% find the stroke plane angle

wing_tip=[squeeze(Q_R(:,2,:)) squeeze(-Q_L(:,2,:))];
h_wing_tip = figure;
plot3(wing_tip(1,:),wing_tip(2,:), wing_tip(3,:),'b.');
set(gca,'YDir','reverse','ZDir','reverse');
xlabel('$\mathbf{b}_x$','interpreter','latex');
ylabel('$\mathbf{b}_y$','interpreter','latex');
zlabel('$\mathbf{b}_z$','interpreter','latex');
axis equal;
grid on;
view(-90,90);

% fit a plane to the points of wing tip
[V D]=eig(wing_tip*wing_tip');
[~, I_min]=min(diag(D));
s=V(:,I_min);
if s(1) < 0
    s=-s;
end
disp(s);
hold on;
plot3([0 s(1)],[0 s(2)],[0 s(3)]);

beta=atan2(-s(3),s(1));
disp(['beta = ' num2str(beta*180/pi)]);

stroke_plane_yaw = atan2(s(2),s(1));
disp(['stroke_plane_yaw_offset = ' num2str(stroke_plane_yaw*180/pi)]);

%% find the wing kinematics angles
% wing attitude is rotated to offset stroke_plane_yaw

for k=1:N
    Stroke_Offset = expm(stroke_plane_yaw*hat(e3));
    Qs_R(:,:,k) = Stroke_Offset'*Q_R(:,:,k);
    Qs_L(:,:,k) = Stroke_Offset'*Q_L(:,:,k);
    [E_R(:,k) E_L(:,k)]= wing_attitude_to_angle(beta, Qs_R(:,:,k), Qs_L(:,:,k));
end

wing_tip=[squeeze(Qs_R(:,2,:)) squeeze(-Qs_L(:,2,:))];
figure;
plot3(wing_tip(1,:),wing_tip(2,:), wing_tip(3,:),'b.');
set(gca,'YDir','reverse','ZDir','reverse');
xlabel('$\mathbf{b}_x$','interpreter','latex');
ylabel('$\mathbf{b}_y$','interpreter','latex');
zlabel('$\mathbf{b}_z$','interpreter','latex');
axis equal;
grid on;
view(-90,90);

% fit a plane to the points of wing tip
[V D]=eig(wing_tip*wing_tip');
[~, I_min]=min(diag(D));
s=V(:,I_min)
if s(1) < 0
    s=-s;
end
hold on;
plot3([0 s(1)],[0 s(2)],[0 s(3)]);

beta=atan2(-s(3),s(1));
disp(['beta = ' num2str(beta*180/pi)]);

stroke_plane_yaw = atan2(s(2),s(1));
disp(['stroke_plane_yaw_offset = ' num2str(stroke_plane_yaw*180/pi)]);


h_E = figure;
my_ylabel={'$\phi$','$\theta$','$\psi$'};
for ii=1:3
    subplot(3,1,ii);
    plot(t,E_R(ii,:)*180/pi,'r');
    hold on;
    plot(t,E_L(ii,:)*180/pi,'b');
    ylabel(my_ylabel{ii},'interpreter','latex');
    plot(t,0.5*(E_R(ii,:)+E_L(ii,:))*180/pi,'k--');    
end
legend('right wing','left wing','average');


%% Fitting with Fourier series

fit_phi=fit(t,180/pi*(E_R(1,:)+E_L(1,:))'/2,'fourier5');
fit_theta=fit(t,180/pi*(E_R(2,:)+E_L(2,:))'/2,'fourier5');
fit_psi=fit(t,180/pi*(E_R(3,:)+E_L(3,:))'/2,'fourier1');
fit_theta_th=fit(t,180/pi*theta_th','fourier1');
fit_theta_ab=fit(t,180/pi*theta_ab','fourier1');

disp('copy to ../wing_kinematics.m');
F_phi=save_fit(fit_phi,5);
F_theta=save_fit(fit_theta,5);
F_psi=save_fit(fit_psi,1);
disp('copy to ../body_attitude.m');
F_theta_th=save_fit(fit_theta_th,1);
disp('copy to ../abdomen_attitude.m');
F_theta_ab=save_fit(fit_theta_ab,1);

filename='fit_VICON_data';
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save(filename,saveVars{:});
evalin('base',['load ' filename]);

bool_print = true;
if bool_print
    figure(h_wing_tip);pause(1);
    print('fit_VICON_wing_tip','-depsc');
    figure(h_E);
    print('fit_VICON_E','-depsc');    
    !mv fit_VICON*.eps ../../doc/Figs/
end

end

function F = save_fit(fit, N_order)
% save the Fourier coefficients into a vector
F.f=fit.w/2/pi;
F.A0=fit.a0;
for i=1:N_order
    eval(['F.AN(' num2str(i) ')=fit.a' num2str(i) ';']);
    eval(['F.BN(' num2str(i) ')=fit.b' num2str(i) ';']);
end

disp(['F.f = ' num2str(F.f) ';']);
disp(['F.A0 = ' num2str(F.A0) ';']);
disp(['F.AN = ' mat2str(F.AN) ';']);
disp(['F.BN = ' mat2str(F.BN) ';']);

end

function [E_R, E_L] = wing_attitude_to_angle(beta, Q_R, Q_L)
e2=[0 1 0]';

exp_beta = expm(beta*hat(e2));

Q_R=exp_beta'*Q_R;
psi_R=asin(Q_R(1,2));
theta_R=atan2(Q_R(1,3)/cos(psi_R), Q_R(1,1)/cos(psi_R));
phi_R=atan2(Q_R(3,2)/cos(psi_R), Q_R(2,2)/cos(psi_R));

Q_L=exp_beta'*Q_L;
psi_L = asin(-Q_L(1,2));
theta_L=atan2(Q_L(1,3)/cos(psi_L), Q_L(1,1)/cos(psi_L));
phi_L=atan2(-Q_L(3,2)/cos(psi_L), Q_L(2,2)/cos(psi_L));

E_R = [phi_R; theta_R; psi_R];
E_L = [phi_L; theta_L; psi_L];
end


