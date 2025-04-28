function [cos_value, L2_error] =  util_similarity(beta_estimated, beta_compare,m)

cos_value = zeros(m,1);
L2_error = zeros(m,1);
for i = 1:m
     
    cos_value(i) = beta_estimated(i,:)*beta_compare(i,:)'/...
        norm(beta_estimated(i,:))/ norm(beta_compare(i,:));
    L2_error(i) = 2-2*abs(cos_value(i));

end

end