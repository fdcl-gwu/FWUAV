function [m, c, Rsq] = lin_reg(x, y)
% Fit a linear model to the given data, return slope and intercept
    X = ones(length(x), 1+size(x,2));
    if size(X, 1) == size(y, 1)
        X(:, 2:end) = x;
        beta = (X'*X) \ (X'*y);
        c = beta(1);
        m = beta(2:end);
        yCalc = X * beta;
        Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
    end
end
