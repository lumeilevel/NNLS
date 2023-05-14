function Q = Q_tau(X, Y, b, tau, mu, Amap, ATmap)
    y = Amap(Y) - b;
    Z = X - Y;
    Q = 0.5*(y'*y) + trace(ATmap(y)'*Z) + tau*0.5*norm(Z, 'fro')^2 + mu*sum(svd(X));
end