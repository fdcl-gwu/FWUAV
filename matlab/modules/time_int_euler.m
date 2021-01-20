% Euler's method:
function [x,fn] = time_int_euler(fun,INSECT,WK_R,WK_L,h,t,x,xDot,x0,SIM)

% Set the initial condition:
x(:,1) = x0;

% One step of Euler's method
for i = 1
    % First call to Stream to get it started:
    xDot(:,i) = fun(INSECT,WK_R,WK_L,t(1),x(:,1),1,SIM);
    fn = xDot(:,i);
    xn = x(:,i);
    x(:,2) = xn + h*fn;
    % This gets us to timeStep = 2.
end
end