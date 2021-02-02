function [U_t, tau_t, theta] = torsion_terms(INSECT, Q, beta)
%% Potential energy from torsional spring
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

sx = Q' * expmhat(beta*e2) * e1;
theta = atan2(e3'*sx, e1'*sx);
U_t = 0.5 * INSECT.Ktorsion * theta^2;
% tau_flex = (INSECT.Ktorsion*theta)/((1+tan(theta)^2) * (e1'*sx)^2) * hat(sx)^2 * e2;
tau_t = (INSECT.Ktorsion*theta)/((1+tan(theta)^2) * (e1'*sx)^2) * ((e2'*sx)*sx - e2);
end
