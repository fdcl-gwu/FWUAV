function [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK)
%wing_kinematics: compute wing Euler angles and their time-derivaties
%
% [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK) computes
%
%       Euler = [phi theta psi]' (flapping, pitch, deviation)
%
%for a given time and a struct variables with the following members
%describing wing kinematics
%
%         WK.f, Wk.beta
%         WK.phi_m, WK.phi_K, WK.phi_0
%         WK.theta_m, WK.theta_C, WK.theta_0, WK.theta_a
%         WK.psi_m, WK.psi_N, WK.psi_a, WK.psi_0
%

switch WK.type
    case 'Monarch'
        
        F_phi.A0 = [6.6842];
        F_phi.AN = [64.5842; 9.4086; -4.3509; -2.2789; 1.8163; 0.7403; 0.0321; 0.0660; 0.2409; 0.2811;];
        F_phi.BN = [-2.9875; 20.7364; 2.5818; -2.7271; 0.0964; 1.0625; 1.0663; 0.5153; 0.4839; 0.5237;];
        
        [phi phi_dot phi_ddot]=Fourier_eval(t, WK.f, F_phi);
        
        F_theta.A0 = [13.7032];
        F_theta.AN = [29.5026; -9.9156; 1.1941; 2.8979; -0.9117; 0.6347; 0.0604; 0.0986; 0.1897; 0.0793;];
        F_theta.BN = [31.1984; 10.5117; -7.3804; 2.9382; -0.1440; 0.3463; 0.4723; 0.2309; 0.3287; 0.2771;];
        
        [theta theta_dot theta_ddot]=Fourier_eval(t, WK.f, F_theta);
        
        F_psi.A0 = [10.9516];
        F_psi.AN = [34.8540; -3.0510; -3.9149; -1.9479; 1.1978; -0.0376; - 0.1545; -0.1020; 0.0494; 0.0893;];
        F_psi.BN = [3.2831; 4.0611; 4.7746; -0.3536; -0.4101; 1.4448; 0.5952; 0.3598; 0.3204; 0.3519;];
        
        
        [psi psi_dot psi_ddot]=Fourier_eval(t, WK.f, F_psi);
%         psi=0;
%         psi_dot=0;
%         psi_ddot=0;
        
    %case 'BermanWang'
    otherwise
        
        %% phi / flapping
        A=WK.phi_m / asin(WK.phi_K);
        a=WK.phi_K;
        b=2*pi*WK.f;
        
        phi = A*asin( a * cos(b*t)) + WK.phi_0;
        phi_dot = -(A*a*b*sin(b*t))/(1 - a^2*cos(b*t)^2)^(1/2);
        phi_ddot = (A*a*b^2*cos(b*t)*(a^2 - 1))/(1 - a^2*cos(b*t)^2)^(3/2);
        
        %% theta / pitching
        A=WK.theta_m / tanh(WK.theta_C);
        a=WK.theta_C;
        b=2*pi*WK.f;
        c=WK.theta_a;
        
        theta = A * tanh( a * sin(b*t + c) ) + WK.theta_0;
        theta_dot = -A*a*b*cos(c + b*t)*(tanh(a*sin(c + b*t))^2 - 1);
        theta_ddot = A*a*b^2*(tanh(a*sin(c + b*t))^2 - 1)*(sin(c + b*t) + 2*a*cos(c + b*t)^2*tanh(a*sin(c + b*t)));
        
        %%  psi / deviation
        A=WK.psi_m;
        a=2*pi*WK.psi_N*WK.f;
        b=WK.psi_a;
        
        psi = A * cos( a*t + b ) + WK.psi_0;
        psi_dot  = A * -a * sin(a*t+b);
        psi_ddot = A * -a^2 * cos(a*t+b);
        
        
end

%% return values
Euler=[phi theta psi]';
Euler_dot=[phi_dot theta_dot psi_dot]';
Euler_ddot=[phi_ddot theta_ddot psi_ddot]';
end

function [a a_dot a_ddot]= Fourier_eval(t, f, F)
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

