function s = vee(S)
% Get the element in R^{3} corresponding to one in so(3)

% if ~all(size(S) == [3, 3])
%     error('Incorrect input to vee map\n');
% end

s=[-S(2,3); S(1,3); -S(1,2)];

end
