function [Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell)
% Y_cell = generate_response(m,n,d,X_cell,e_cell,xi_cell)
% Generate the response variable according to the sample and noise
% m - number of client
% n - sample size (m×1)
% d - dimension of samples
% real_beta - real coefficents on each client (m×d)
% X_cell - samples with (n×d)
% e_cell - noise (n×1)
% xi_cell - sign flip variable (n×1)

Y_cell = cell(m,1);
YN = zeros(sum(n),1);

ni=1;
for i = 1:m
    
    nn = n(i);
    e = cell2mat(e_cell(i,1));
    xi = cell2mat(xi_cell(i,1));
    X = cell2mat(X_cell(i,1));
    beta = real_beta(i,:)'; %(d×1)
    
    Y = xi.*sign(X*beta+e);
    
    Y_cell(i,1) = mat2cell(Y,nn,1);
    
    YN(ni:ni+n(i)-1) = Y;
    ni = n(i)+1;
end




end