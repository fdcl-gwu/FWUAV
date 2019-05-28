function fit_exp_data
%FIT_EXP_DATA Load experimental data and fit into Fourier series
%TODO: least-square fit for the same frequency consistent for all of wing
%kinematics angles and pitch angles.:
clear all;
close all;

XLS=readtable('Time_history_data_for_Dr.Lee.xlsx');
XLS=table2array(XLS);
%%
t=XLS(:,1)-XLS(1,1);
phi = XLS(:,15);
theta = -XLS(:,16); % the convention of the wing pitch angle is opposite 
psi = XLS(:,17);

fit_phi=fit(t,phi,'fourier5');
fit_theta=fit(t,theta,'fourier5');
fit_psi=fit(t,psi,'fourier5');

%% the first cycle of theta_th and theta_ab are irregular and excluded for fitting
I_fit=21:length(XLS);

t=XLS(I_fit,1)-XLS(1,1);
theta_th=XLS(I_fit,8);
theta_ab=XLS(I_fit,9)-XLS(I_fit,8); % the abdoment pitch angle, namely theta_ab is defined as the angle RELATIVE to the thorax

fit_theta_th=fit(t,theta_th,'fourier5');
fit_theta_ab=fit(t,theta_ab,'fourier1');

F_phi=save_fit(fit_phi,5);
F_theta=save_fit(fit_theta,5);
F_psi=save_fit(fit_psi,5);
F_theta_th=save_fit(fit_theta_th,5);
F_theta_ab=save_fit(fit_theta_ab,1);

filename='fit_exp_data';
save(filename);
evalin('base',['load ' filename]);

end

function F = save_fit(fit, N_order)
% save the Fourier coefficients into a vector
F.f=fit.w/2/pi;
F.A0=fit.a0;
for i=1:N_order
    eval(['F.AN(' num2str(i) ')=fit.a' num2str(i) ';']);
    eval(['F.BN(' num2str(i) ')=fit.b' num2str(i) ';']);
end

disp(F);
end

