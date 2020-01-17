function [a, a_dot, a_ddot]= eval_Fourier(t, f, F)
% Evaluate the function with given fourier coefficients
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
