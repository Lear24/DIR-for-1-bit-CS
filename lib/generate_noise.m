function [e_cell, xi_cell] = generate_noise(m,n,sigma,q)
% [e_cell, xi_cell] = generate_noise(m,n,sigma,q)
% Assign noise on each client according to sample size and noise variance
% m - number of clients
% n - m-varibles (m×1)
% sigma, q - the variance of noise and the probability of sign flipping in each client

e_cell = cell(m,1); xi_cell = cell(m,1);
for i =1:m
    nn = n(i);
    sig = sigma(i);
    qq = q(i);
    e = randn(nn,1)*sig;
    xi=2*binornd(1,qq,nn,1)-1;
    
    e_cell(i,1)=mat2cell(e,nn,1);
    xi_cell(i,1)=mat2cell(xi,nn,1);
end




end