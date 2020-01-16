function [R W W_dot theta]=body_attitude(t, f, WK)
%body_attitude compute attitude of body
%
% [R W W_dot theta]=body_attitude(t, f) returns the attitude, the angular
% velocity, the angular acceleration, and the pitch angle of body, computed
% according to the fitted experimental data (see exp_data/fit_exp_data.m)
%
%   INPUT
%       t:  time (sec)
%       f:  flapping frequench (Hz)
%   OUTPUT
%       R:  attitude of the body
%       W:  angular velocity of the body (rad/sec)
%       W_dot:  angular acceleration of the body (rad/sec^2)
%       theta:  pitch angle of the body (rad)
%
% [R W W_dot]=body_attitude(theta) returns a fixed attitude for the given
% body pitch angle
%
%   INPUT
%       theta:  pitch angle of the body (rad)
%   OUTPUT
%       R:  attitude of the body
%       W:  angular velocity of the body (rad/sec)
%       W_dot:  angular acceleration of the body (rad/sec^2)

e2=[0 1 0]';

switch WK.bo_type
    case 'fixed'
        % fixed body attidue
        theta=WK.theta_B_0;
        R=expm(theta*hat(e2));
        W=zeros(3,1);
        W_dot=zeros(3,1);    
    case 'experimental'
%         t=varargin{1};
%         f=varargin{2};
%         Data constructed by ./exp_data/fit_exp_data.m    
%         F_theta_th.f=10.2468;
%         F_theta_th.A0=19.7703;
%         F_theta_th.AN=[0.7632 2.1236 -0.1915 -0.3702 -0.1729];
%         F_theta_th.BN=[9.3219 0.4207 -0.9178 -0.7781 0.3443];
% 
%         [theta theta_dot theta_ddot]= eval_Fourier(t, f, F_theta_th);
% 
%         Data constructed by ./exp_data/fit_VICON_data.m    
        F_theta_th.f = 10.2213;
        F_theta_th.A0 = 18.6094-3;
        F_theta_th.AN = 1.1952126587572*1.3;
        F_theta_th.BN = 8.21314837914776*1.3;
        [theta theta_dot theta_ddot]= eval_Fourier(t, f, F_theta_th);

        R=expm(theta*hat(e2));
        W=theta_dot*e2;
        W_dot=theta_ddot*e2;
    case 'varying'
%         t=varargin{1};
%         WK=varargin{2};
        A=WK.theta_B_m;
        a=2*pi*WK.f;
        b=WK.theta_B_a;
        theta = A * cos( a*t + b ) + WK.theta_B_0;
        theta_dot  = A * -a * sin(a*t+b);
        theta_ddot = A * -a^2 * cos(a*t+b);  
        R=expm(theta*hat(e2));
        W=theta_dot*e2;
        W_dot=theta_ddot*e2;
end

end

function [a a_dot a_ddot]= eval_Fourier(t, f, F)
a=F.A0;
a_dot=0;
a_ddot=0;

for q = 1:length(F.AN)
    a = a + F.AN(q)*cos(2*pi*q*f*t) + F.BN(q)*sin(2*pi*q*f*t);
    a_dot = a_dot -(2*pi*q*f)*F.AN(q)*sin(2*pi*q*f*t) + ...
        (2*pi*q*f)*F.BN(q)*cos(2*pi*q*f*t);
    a_ddot = a_ddot -(2*pi*q*f)^2*F.AN(q)*cos(2*pi*q*f*t) - ...
        (2*pi*q*f)^2*F.BN(q)*sin(2*pi*q*f*t);    
end

a=a*pi/180;
a_dot=a_dot*pi/180;
a_ddot=a_ddot*pi/180;
end
