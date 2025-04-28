function [sigma,q] = assign_noise_variance(m,sig,qq,ratio)
% [sigma,q] = assign_noise_variance(m,sig,qq,ratio2) 
% Assign the amount of noise on each client
% m - number of client
% sig - two types of noise variance, [sigma1,sigma2]
% qq - two types of flipping probability, [q1,q2]
% ratio - the ratio of clients with [sigma1,q1]
sigma =zeros(m,1);q =zeros(m,1);
m1 = round(ratio*m);
mm = zeros(m,1);
mm(1:m1,1) = ones(m1,1);
index = randperm(m,m);
mm=mm(index); %随机分配的噪声
for i = 1:m
    if mm(i) == 1
        
        sigma(i) = sig(1);
        q(i) = qq(1);
        
    elseif mm(i) == 0
        
        sigma(i) = sig(2);
        q(i) = qq(2);
        
    end
end
end