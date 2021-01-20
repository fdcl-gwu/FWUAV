% 2nd order Adams-Bashforth method:
function x = time_int_adam_bash(fun,INSECT,WK_R,WK_L,t,x,xDot,x0,SIM)
h = t(2)-t(1);
% Get first two time steps:
% First is from initial conditions
% Second is from one application of Euler's method
[x(:,1:2), xDot(:,1)] = time_int_euler(fun,INSECT,WK_R,WK_L,h,t(1),x(:,1),xDot(:,1),x0,SIM);

% Now to integrate from timestep = 3:end
for i = 1:length(t)-1
    fn = xDot(:,i);
    xnp1 = x(:,i+1);
    xDot(:,i+1) = fun(INSECT,WK_R,WK_L,t(i+1),x(:,i+1),i+1,SIM);
    fnp1 = xDot(:,i+1);
    if i <= length(t)-2
        x(:,i+2) = xnp1 + h*(3/2*fnp1 - 1/2*fn);
    end
end

end