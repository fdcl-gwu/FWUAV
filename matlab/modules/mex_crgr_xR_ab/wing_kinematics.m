function [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK)
coder.inline('always');
%wing_kinematics: compute wing kinematics angles and their time-derivaties
%
% [Euler, Euler_dot, Euler_ddot] = wing_kinematics(t,WK) computes the wing
% kinematics angles defined by
%
%       Euler = [phi theta psi]' (flapping, pitch, deviation)
%
% and their derivatives for a given time t, and the struct variables WK
% with the following members 
%
%         WK.f     flapping frequency
%         Wk.beta  stroke plane angle
%         WK.type  wing kinematics type
%                   "Monarch"       wing kinematics from Monarch experiments
%                   "BermanWang"    a min energy model proposed by Berman
%                   and Wang with the following variables (see Section 1.3)
%
%         WK.phi_m, WK.phi_K, WK.phi_0
%         WK.theta_m, WK.theta_C, WK.theta_0, WK.theta_a
%         WK.psi_m, WK.psi_N, WK.psi_a, WK.psi_0
%

% switch WK.type
%     case 'Monarch'
%                
%         % Data constructed by ./exp_data/fit_VICON_data.m
%         F_phi.f = 10.2247;
%         F_phi.A0 = -0.83724;
%         F_phi.AN = [59.7558745992612 -0.137466978762473 0.137978226025185 0.433746159939459 -0.074830919096204];
%         F_phi.BN = [-6.55441083419535 6.435543953825 -2.05072909120033 0.221239708063663 0.0575790280561444];
%         F_theta.f = 10.1838;
%         F_theta.A0 = -5.9583;
%         F_theta.AN = [24.0863053541935 -4.46860932682531 8.04218451917262 -0.601926108817941 -0.642171907559121];
%         F_theta.BN = [30.5848624271628 4.64142325712464 -3.83504323967398 2.01570854870463 -1.70930433046515];
%         F_psi.f = 20.2767;
%         F_psi.A0 = -1.9753;
%         F_psi.AN = -2.063160749463;
%         F_psi.BN = -2.60158935092774;
%         [phi phi_dot phi_ddot]=eval_Fourier(t, WK.f, F_phi);
%         [theta theta_dot theta_ddot]=eval_Fourier(t, WK.f, F_theta);
%         [psi psi_dot psi_ddot]=eval_Fourier(t, 2*WK.f, F_psi);
% 
%         %case 'BermanWang'
%     otherwise
        % phi / flapping
        A=WK.phi_m / asin(WK.phi_K);
        a=WK.phi_K;
        b=2*pi*WK.f;
        
        phi = A*asin( a * cos(b*t)) + WK.phi_0;
        phi_dot = -(A*a*b*sin(b*t))/(1 - a^2*cos(b*t)^2)^(1/2);
        phi_ddot = (A*a*b^2*cos(b*t)*(a^2 - 1))/(1 - a^2*cos(b*t)^2)^(3/2);
        
        % theta / pitching
        A=WK.theta_m / tanh(WK.theta_C);
        a=WK.theta_C;
        b=2*pi*WK.f;
        c=WK.theta_a;
        
        theta = A * tanh( a * sin(b*t + c) ) + WK.theta_0;
        theta_dot = -A*a*b*cos(c + b*t)*(tanh(a*sin(c + b*t))^2 - 1);
        theta_ddot = A*a*b^2*(tanh(a*sin(c + b*t))^2 - 1)*(sin(c + b*t) + 2*a*cos(c + b*t)^2*tanh(a*sin(c + b*t)));
        
        %  psi / deviation
        A=WK.psi_m;
        a=2*pi*WK.psi_N*WK.f;
        b=WK.psi_a;
        
        psi = A * cos( a*t + b ) + WK.psi_0;
        psi_dot  = A * -a * sin(a*t+b);
        psi_ddot = A * -a^2 * cos(a*t+b);               
% end

%% return values
Euler=[phi theta psi]';
Euler_dot=[phi_dot theta_dot psi_dot]';
Euler_ddot=[phi_ddot theta_ddot psi_ddot]';
end
