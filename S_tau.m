function [SG, sd] = S_tau(G, tau, mu)
    [U,S,V] = svd(G, 'econ');
%     [U, S, V] = mexsvd(G);
    sd = max(S - mu/tau, 0);
    SG = U * sd * V';
    sd = diag(sd);
end