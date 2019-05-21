function [R W W_dot]=body_attitude(t,varargin)
e2=[0 1 0]';

if nargin >=2
    bool_fixed=varargin{1};
end

if ~bool_fixed
    % compute attitude as a function of time here
else
    % fixed body attidue 
    R=expm(25*pi/180*hat(e2));
    W=zeros(3,1);
    W_dot=zeros(3,1);
end

end
