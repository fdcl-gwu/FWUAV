function [Q_A W_A W_A_dot theta_A]=abdomen_attitude(t, f, WK)
%abdomen_attitude compute attitude of abdomen relative to the body
%
% [Q_A W_A W_A_dot theta_A]=abdomen_attitude(t, f) returns the attitude, the angular
% velocity, the angular acceleration, and the pitch angle of abdomen, computed
% according to the fitted experimental data (see exp_data/fit_exp_data.m)
%
%   INPUT
%       t:  time (sec)
%       f:  flapping frequench (Hz)
%   OUTPUT
%       Q_A:  attitude of the abdomen relative to the body
%       W_A:  relative angular velocity of the abdomen (rad/sec)
%       W_A_dot:  relative acceleration of the abdomen (rad/sec^2)
%       theta_A:  relative pitch angle of the abdomen (rad)
%
% [Q_A W_A W_A_dot]=abdomen_attitude(theta_A) returns a fixed attitude for the
% given body pitch angle
%
%   INPUT
%       theta_A:  relative pitch angle of the abdomen (rad)
%   OUTPUT
%       Q_A:  attitude of the abdomen relative to the body
%       W_A:  relative angular velocity of the abdomen (rad/sec)
%       W_A_dot:  relative acceleration of the abdomen (rad/sec^2)

e2=[0 1 0]';

switch WK.ab_type
    case 'fixed'
        % fixed body attidue
        theta_A=WK.theta_A_0;
        Q_A=expm(theta_A*hat(e2));
        W_A=zeros(3,1);
        W_A_dot=zeros(3,1);
    case 'experimental'
%         t=varargin{1};
%         f=varargin{2};
        % Data constructed by ./exp_data/fit_VICON_data.m
        F_theta_ab.f = 10.1756;
        F_theta_ab.A0 = 17.3417;
        F_theta_ab.AN = -13.2134053462198;
        F_theta_ab.BN = -3.91351009294982;

        [theta_A theta_A_dot theta_A_ddot]= eval_Fourier(t, f, F_theta_ab);

    %     % abdomen with the opposite phase
    %     theta_A = 0.6053 - theta_A;
    %     theta_A_dot = -theta_A_dot;
    %     theta_A_ddot = -theta_A_ddot;

        Q_A=expm(theta_A*hat(e2));
        W_A=theta_A_dot*e2;
        W_A_dot=theta_A_ddot*e2;
    case 'varying'
%         t=varargin{1};
%         WK=varargin{2};
        A=WK.theta_A_m;
        a=2*pi*WK.f;
        b=WK.theta_A_a;
        theta_A = A * cos( a*t + b ) + WK.theta_A_0;
        theta_A_dot  = A * -a * sin(a*t+b);
        theta_A_ddot = A * -a^2 * cos(a*t+b);  
        Q_A=expm(theta_A*hat(e2));
        W_A=theta_A_dot*e2;
        W_A_dot=theta_A_ddot*e2;
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
