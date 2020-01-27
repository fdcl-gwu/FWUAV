function R = expmhat(W)
% Exponential of an element in Lie algebra of SO(3) dentified by R^{3}
% using Rodrigues formula.
    theta = norm(W);
    if theta == 0
        R = eye(3);
    else
        W = W / theta;
        R = eye(3) + sin(theta) * hat(W) + (1 - cos(theta)) * (hat(W)^2);
    end
end
