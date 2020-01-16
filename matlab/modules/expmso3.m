function A=expmso3(s)
% Compute the element in SO(3) corresponding to an element in its Lie
% algebra identified by R^{3}
A = expm(hat(s));
end
