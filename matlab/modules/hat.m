function hat_x = hat(x)
% Get the element in so(3) corresponding to one in R^{3}

% if numel(x) ~= 3
%     error('Incorrect input to hat map\n');
% end

hat_x=[0 -x(3) x(2);
    x(3) 0 -x(1);
    -x(2) x(1) 0];

end
