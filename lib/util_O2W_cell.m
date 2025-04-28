function [W_cell,norm_W_vec]=util_O2W_cell(beta_cell,m,d)
% [W_cell,norm_W_vec]=util_O2W_cell(beta_cell,m,d)
% 
% On all clinets. tranfer the coefficents in problem (O) to the corresponding coefficents in problem (W).
% beta_cell: coefficient in (O) in form of (m×1) cell.
% 
% W: the  (d+1 × d+1) variable define in (W).
% norm_W: the norm of beta, 
% norm_W_1: trace(W'W)-1, equivalently the norm of beta when W is rank-1 matrix.

W_cell = cell(m,1);
norm_W_vec = zeros(m,1);

for i=1:m
    beta=cell2mat(beta_cell(i));
    W = [beta';1]*[beta,1];
    W_cell(i)=mat2cell(W,d+1);
    
    
    % 此处的norm_W和tr(A'A)-1 (norm_W_1) 一定是相等的.
    norm_W = norm(beta,2)^2;
    norm_W_vec(i) = norm_W;
    
    
    % 检验 其实放在这里大概率没用, 可以之后再用上.
    [V,D] = eig(W);
    D(D > -0.0001 & D < 0)=0;
    D(D < 0.0001 & D > 0)=0;
    L=V*(D^0.5);
    norm_W_1 = trace(L'*L)-1;
    if abs(norm_W-norm_W_1)>10^(-4)
        disp('norm计算出现错误, 两种计算方式不匹配')
        disp(norm_W-norm_W_1)
        pause;
    end
end


end