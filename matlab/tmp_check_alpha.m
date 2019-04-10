function tmp_check_alpha
close all;

N=200;
theta=linspace(0,2*pi,N);

eps=1e-6;

norm_U=1.531234;
for k=1:N
    U(:,k)=[cos(theta(k)), 0, sin(theta(k))]*norm_U;
    delU(:,k)=[rand 0 rand]';
    [alpha(k) alpha_dot(k)]=compute_alpha(U(:,k),delU(:,k));
    
    alpha_new(k)=compute_alpha(U(:,k)+eps*delU(:,k));
    
end

plot((alpha_new-alpha)/eps,'r');
hold on;
plot(alpha_dot,'b--');


filename='tmp_check_alpha';
save(filename);
evalin('base',['load ' filename]);


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
    
    if alpha < pi/4
        alpha_dot = -abs(U(3))*dot(U,U_dot)/norm_U^3 + sign(U(3))*U_dot(3)/norm_U;
        alpha_dot = alpha_dot / cos(alpha);       
    else
        alpha_dot = abs(U(1))*dot(U,U_dot)/norm_U^3 - sign(U(1))*U_dot(1)/norm_U;
        alpha_dot = alpha_dot / sin(alpha);
    end
       
    varargout{1}=alpha_dot;
end

end
