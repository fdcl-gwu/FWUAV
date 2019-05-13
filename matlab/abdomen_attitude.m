function [Q_A W_A W_A_dot]=abdoment_attitude(t,varargin)

if nargin >=2
    bool_fixed=varargin{1};
end

if ~bool_fixed
    % compute abdoment attitude as a function of time here
else
    Q_A=eye(3);
    W_A=zeros(3,1);
    W_A_dot=zeros(3,1);
end

end


