function value = f(X, b, mu, Amap)
    x = Amap(X) - b;
    value = 0.5 * (x' * x) + mu * sum(svd(X));
% cvx_begin
%     value = 0.5 * (x' * x) + mu * norm_nuc(X);
% cvx_end
end