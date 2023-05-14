function [SG, sd] = S_tau_truncation(G, tau, mu, gap)
    [U,S,V] = svd(G, 'econ');
    rho = diag(max(S - mu/tau, 0));
    q = sum(rho > 0);
    for j = 1 : q - 1
        chi = mean(rho(1:j)) / mean(rho(j+1:q));
        if chi >= gap
            rho(j+1:q) = 0;
            break
        end
    end
    sd = diag(rho);
    SG = U * sd * V';
end