clear;
n = 1000;
kk = [5,10,15,20,25,30,35,40,45,50];
kk = kk*10;
% A = sprandsym(n, n);
A = randn(n, n);
A = A + A';
L = length(kk);
time1 = zeros(L, 1);
time2 = zeros(L, 1);
time3 = zeros(L, 1);
for i = 1 : L
    k = kk(i);
    tic;
    [U1, S1, V1] = lansvd(A, k, 'L');
    time1(i) = toc;

    tic;
    [V, D] = eigs(A*A', k);
    [U2, S2, V2] = svd(V'*A, 'econ');
    U2 = V * U2;
    time2(i) = toc;

    tic;
    [U3, S3, V3] = svd(A);
    time3(i) = toc;
end
plot(kk,time1,'-+',kk,time2,'-*',kk,time3,'-x');
legend('lansvd','approx eig','fullsvd','Location','Northwest');
title('The cpu time for different methods');
xlabel('The predetermined rank');
ylabel('cpu time /s');