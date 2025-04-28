function beta_ls = onbitCS_LS(X_cell,Y_cell,m,n,d)


beta_ls = zeros(d,m);

for i = 1:m
    
    X = cell2mat(X_cell(i));
    Y = cell2mat(Y_cell(i));
    beta_ls(:,i) = (X'*X)^(-1)*X'*Y;
    
end


beta_ls = beta_ls';



end