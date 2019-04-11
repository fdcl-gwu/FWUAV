function [alpha alpha_dot U_R U_theta]=wing_QS_aerodynamics(MONARCH, W_R, W_R_dot)
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

U_R = MONARCH.tilde_r_2*MONARCH.l*cross(W_R,e2);
U_R_dot = MONARCH.tilde_r_2*MONARCH.l*cross(W_R_dot,e2);
[alpha alpha_dot U_theta]=compute_alpha(U_R, U_R_dot);

end

function [alpha varargout]=compute_alpha(U, varargin)
% compute the angle of attack
% input
%       U : velocity of the wind projected to the x-z plane of the wing
%       frame
%       U_dot : (optional) time-derivative of U
% output
%       alpha : angle of attack [0,pi/2]
%       alpha_dot : computed when U_dot is provided

norm_U=norm(U);
alpha=acos(abs(U(1))/norm_U);

if nargin > 1
    tmp=zeros(1,2);
    U_dot=varargin{1};
    
    %if alpha < pi/4
        alpha_dot(1) = -abs(U(3))*dot(U,U_dot)/norm_U^3 + sign(U(3))*U_dot(3)/norm_U;
        alpha_dot(1) = alpha_dot(1) / cos(alpha);       
    %else
        alpha_dot(2) = abs(U(1))*dot(U,U_dot)/norm_U^3 - sign(U(1))*U_dot(1)/norm_U;
        alpha_dot(2) = alpha_dot(2) / sin(alpha);
    %end
    if norm(alpha_dot(1)-alpha_dot(2)) > 10
        keyboard;
    end
    
    varargout{1}=alpha_dot;
end

varargout{2}=atan2(U(3),U(1));

end

