function R = expmhat(W)
% Exponential of an element in Lie algebra of SO(3) dentified by R^{3}
% using Rodrigues formula. This is faster than using {expm}.
    theta = norm(W);
    if theta ~= 0
        W = W / theta;
        hat_W = hat(W);
        R = eye(3) + sin(theta) * hat_W + (1 - cos(theta)) * (hat_W * hat_W);
    else
        R = eye(3);
    end
end
