function [X, i, time, hist] = myAPGL1(nr, nc, Amap, ATmap, b, mu_target, Lf, eta, tol)
    tic;
    X0 = zeros(nr, nc);
    X = X0;
    t0 = 1;
    t = t0;
    tau = Lf;
    max_iter = 200;
    normb = max(1, norm(b));
    objold = normb * mu_target;
    mu = mu_target * 1e4;
%     mu = mu_target;
    for i = 1 : max_iter
        Y = X + (t0-1)/t * (X-X0);
        
        tau = eta * tau;
%         for j = epoch
        while 1
            G = Y - 1 / tau * ATmap(Amap(Y)-b);
            S_tau_G = S_tau(G, tau, mu);
            if f(S_tau_G, b, mu, Amap) <= Q_tau(S_tau_G, Y, b, tau, mu, Amap, ATmap)
                break
            else 
                tau = min(Lf, tau / eta);
            end
        end
        
        X0 = X;
        X = S_tau(G, tau, mu);
%         X = S_tau_truncation(G, tau, mu, 100);
        t0 = t;
        t = 0.5*(1 + sqrt(1+4*t*t));
        
%         S = tau*(Y-X)+ATmap(Amap(X-Y));
        obj = f(X, b, mu_target, Amap);
        hist.obj(i) = obj;
        mu = max(mu_target, 0.7*mu);
        hist.relRes(i) = norm(Amap(X)-b) / normb;
        hist.relDist(i) = norm(X-X0,'fro')/norm(X,'fro');
        hist.relObjdiff(i) = abs(obj - objold) / max(obj, 1);
        objold = obj;
%         if norm(S, 'fro') / (tau*max(1, norm(X, 'fro'))) <= tol
        if (mu == mu_target) && (i > 4)
            if max(hist.relDist(i), 0.1*hist.relObjdiff(i)) < tol
                fprintf("\n APGL Early Stopping--iteration: %d\n", i);
                fprintf('[a] relDist < %3.2e', tol);
                fprintf("norm(X-Xold,'fro')/norm(X,'fro') = %f\n", hist.relDist(i));
                break
            end
            if max(0.5*hist.relDist(i), 100*hist.relObjdiff(i)) < tol
                fprintf("\n APGL Early Stopping--iteration: %d\n", i);
                fprintf('[b] relObjdiff < %3.2e', 0.01*tol);
                break;
            end
            idx = i-2:i;
            relResratio = hist.relRes(idx-1)./ hist.relRes(idx);
            if (all(abs(relResratio-1) < 5*tol) || (hist.relObjdiff(i) < 0.1*tol)) && (hist.relDist(i) < 10*tol)
                msg = sprintf('[c] lack of progress in relRes');
                fprintf('\n %s : %2.1e',msg,max(relResratio-1)); 
                break; 
            end
        end
    end
    time = toc;
end