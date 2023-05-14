function [X, i, time, hist] = myAPGL3(nr, nc, Amap, ATmap, b, mu_target, Lf, eta, tol)
    tic;
    X0 = zeros(nr, nc);
    X = X0;
    t0 = 1;
    t = t0;
    tau = Lf;
    max_iter = 1000;
%     mu = mu_target * 1e4;
    mu = mu_target;
    for i = 1 : max_iter
        Y = X + (t0-1)/t * (X-X0);
        G = Y - 1 / tau * ATmap(Amap(Y)-b);
        
%         tau = eta * tau;
%         for j = epoch
%         while 1
%             G = Y - 1 / tau * ATmap(Amap(Y)-b);
%             S_tau_G = S_tau(G, tau, mu);
%             if f(S_tau_G, b, mu, Amap) <= Q_tau(S_tau_G, Y, b, tau, mu, Amap, ATmap)
%                 break
%             else 
%                 tau = min(Lf, tau / eta);
%             end
%         end
        
        X0 = X;
%         X = S_tau(G, tau, mu);
        X = S_tau_truncation(G, tau, mu, 475);
        t0 = t;
        t = 0.5*(1 + sqrt(1+4*t*t));
        
        S = tau*(Y-X)+ATmap(Amap(X-Y));
        hist.obj(i) = f(X, b, mu_target, Amap);
%         mu = max(mu_target, 0.7*mu);
        if norm(S, 'fro') / (tau*max(1, norm(X, 'fro'))) <= tol
            fprintf("\n APGL Early Stopping--iteration: %d\n", i);
            break
        end
    end
    time = toc;
end