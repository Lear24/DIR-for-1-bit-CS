function [X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv)
% X = generation_sample(m,n,d,nu)
% m - number of client
% n - sample size (m×1)
% d - dimension of samples
% nv - parameter for covariance matrix of X
% 输出cell 形式的样本
X_cell = cell(m,1);

% 生成指定的协方差矩阵
Sigma_X = zeros(d,d);
for i = 1:d
    for j = 1:d        
        Sigma_X(i,j) = nv^abs(i-j);  
    end  
end

XN = zeros(sum(n),d);
nn = 1;
for i =1:m
    
   X = mvnrnd(zeros(d,1),Sigma_X,n(i));
   X_cell(i,1) = mat2cell(X,n(i),d);
   XN(nn:nn+n(i)-1,:) = X;
   nn = n(i) + 1;
end

end