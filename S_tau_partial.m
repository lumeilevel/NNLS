function [SG, sd] = S_tau_partial(G, tau, mu, sv)
    [U, S, V] = lansvd(G, sv, 'L');
    
%     [U,S,V] = svd(G, 'econ');

%     [Q, ~] = eigs(G*G', sv);
%     [U, S, V] = svd(Q'*G, 'econ');
%     U = Q * U;
    
    sd = max(S - mu/tau, 0);
    SG = U * sd * V';
    sd = diag(sd);
end