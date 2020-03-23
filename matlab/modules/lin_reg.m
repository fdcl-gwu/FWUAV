function [m, c] = lin_reg(x, y)
% Fit a linear model to the given data, return slope and intercept
    X = ones(length(x), 1+size(x,2));
    if size(X, 1) == size(y, 1)
        X(:, 2:end) = x;
        m = (X'*X) \ (X'*y);
        c = m(1);
        m = m(2:end);
    end
end
