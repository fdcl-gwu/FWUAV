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

        F_phi.A0 = [ 0.8691];
        F_phi.AN = [ 57.8139 5.9774 -1.9430 0.7033 -0.0070 0.2849 0.1503 0.1737 0.2031 0.1588 ];
        F_phi.BN = [ -8.4931 12.5083 -0.4286 0.8218 1.0618 0.3111 0.8475 0.4781 0.5475 0.4647 ];
        
        [phi phi_dot phi_ddot]=Fourier_eval(t+WK.t_shift, WK.f, F_phi);
        
        F_theta.A0 = [11.8841 ];
        F_theta.AN = [-28.1511 5.9864 -1.6381 -1.8261 1.1729 -0.7694 0.0003 0.0416 -0.0818 -0.0119 ];
        F_theta.BN = [-32.2038 -10.1697 5.2429 -3.5196 0.8692 -0.3601 -0.5034 -0.2680 -0.3013 -0.2787 ];
        
        [theta theta_dot theta_ddot]=Fourier_eval(t+WK.t_shift, WK.f, F_theta);
        theta = - theta;
        theta_dot = - theta_dot;
        theta_ddot = - theta_ddot;
        
        F_psi.A0 = [ 16.0238 ];
        F_psi.AN = [ 17.5325 -6.2925 0.3741 -0.5870 0.2648 0.2025 0.2335 0.0750 0.0578 0.1054 ];
        F_psi.BN = [ -5.4048 -1.3062 2.3158 0.9966 0.3479 0.4129 0.1933 0.3538 0.1872 0.2095 ];
                
        [psi psi_dot psi_ddot]=Fourier_eval(t+WK.t_shift, WK.f, F_psi);
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

for q = 1:length(F.AN)-4
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

