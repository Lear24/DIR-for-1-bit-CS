function [beta_E,cos_value,cos_mean]=...
    onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max)
% The distributed decoding Lasso with \lam = 0
% X_cell: (m×1)cell with (n×d)
% Y_cell: (m×1)cell with (n×1)
%
% m: number of clinet
% n: (m×1) vector. The sample size on each clients is different.
% d: dimension.
%
% T_max: The max iteration of conmunication round.
%


%initial-------------------------------------------------------------------
for i = 1
    
    X = cell2mat(X_cell(i));
    Y = cell2mat(Y_cell(i));
    beta_E = (X'*X)^(-1)*X'*Y;
    
    Xbar = X - mean(X);
    bSig1 = Xbar'*Xbar/n(i);
    
end

bSig = zeros(d,d);
Zn = zeros(d,1);
for i = 1:m
    
    nn = n(i);
    X = cell2mat(X_cell(i));
    Y = cell2mat(Y_cell(i));
    
    Xbar = X - mean(X);
    Ybar = Y - mean(Y);
    
    bSig = bSig + Xbar'*Xbar/nn;
    Zn = Zn + Xbar'*Ybar/nn;
    
end
bSig = bSig/m; Zn = Zn/m;
cos_mean = zeros(T_max,1);
cos_value = zeros(m,1);

T = 2;


while T <= T_max
    b = Zn + (bSig1 - bSig)*beta_E;
    beta_E = (bSig1^(-1))*b;
    
    for i = 1:m
        
        cos_value(i) = real_beta(i,:)*beta_E/norm(beta_E)/norm(real_beta(i,:));
        
    end
    
    cos_mean(T) = mean(abs(cos_value));
    
    if T>1 && abs(cos_mean(T) - cos_mean(T-1)) < 10^(-7)
        fprintf('在轮次%d时, 估计效果cos value为 %f \n',T,cos_mean(T));
        cos_mean(T+1:end) = cos_mean(T);
        break;
    end
    if T == T_max
        fprintf('在轮次T_max = %d 未达到收敛, 增加估计轮次500轮 \n',T_max);
        T_max = T_max+500;
    end
    T = T + 1;
end






% for T = 2:T_max
%     b = Zn + (bSig1 - bSig)*beta_E;
%     beta_E = (bSig1^(-1))*b;
%     
%     for i = 1:m
%         
%         cos_value(i) = real_beta(i,:)*beta_E/norm(beta_E)/norm(real_beta(i,:));
%         
%     end
%     
%     cos_mean(T) = mean(abs(cos_value));
%     
%     if T>10 && abs(cos_mean(T) - cos_mean(T-10)) < 10^(-7)
%         fprintf('在轮次%d时, 估计效果cos value为 %f \n',T,cos_mean(T));
%         cos_mean(T+1:end) = cos_mean(T);
%         break;
%     end
%     if T == T_max
%         fprintf('在轮次T_max = %d 未达到收敛, 增加估计轮次500轮 \n',T_max);
%         T_max = T_max+500;
%     end
% end






end



























