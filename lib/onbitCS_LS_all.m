function [beta_ls_all,cos_value1] = onbitCS_LS_all(XN,YN,real_beta,m)


beta_ls_all = (XN'*XN)^(-1)*XN'*YN;
cos_value1 = zeros(m,1);
for i = 1:m
    
    cos_value1(i) = real_beta(i,:)*beta_ls_all/norm(beta_ls_all)/norm(real_beta(i,:));
    
end


end