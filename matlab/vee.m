function s = vee(S)

% if ~all(size(S) == [3, 3])
%     error('Incorrect input to vee map\n');
% end

s=[-S(2,3); S(1,3); -S(1,2)];

end