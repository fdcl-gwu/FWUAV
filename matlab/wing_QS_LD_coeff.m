function [C_L C_D] = wing_QS_LD_coeff(alpha)
%wing_QS_LD_coeff: compute C_L and C_D
%[C_L C_D] = trans_force_coeff(alpha) computes the lift coefficient and the
%drag coefficient for a given angle of attack, alpha in RADIAN

% S. P. Sane and M. H. Dickinson, "The aerodynamic effects of wing rotation
% and a revised quasi-steady model of flapping flight," 
% Journal of experimental biology, vol. 205, no. 8, pp. 1087?1096, 2002.

% convert radian to degree
alpha_deg=alpha*180/pi;

C_L = 0.225 + 1.58 * sind(2.13*alpha_deg -7.2);
C_D = 1.92 - 1.55 * cosd(2.04*alpha_deg-9.82);
end
