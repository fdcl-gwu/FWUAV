function [L_R D_R M_R F_rot_R M_rot_R]=wing_QS_aerodynamics(MONARCH, W_R)
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

U_R = MONARCH.tilde_r_2*MONARCH.l*cross(W_R,e_2);





end

function alpha = angle_of_attack(U)



end
