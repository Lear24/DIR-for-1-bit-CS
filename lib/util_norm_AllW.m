function  norm_AllW = util_norm_AllW(W_cell,m,d)
% Calculates the norm for each parameter for normalization.
% 
% norm_AllW: (m×1) , record the norm of all parameters in (W).
% 
% W_cell: (m×1)cell with (d+1 × d+1) matrix
% m: number of client
% d: dimension of coefficient



norm_AllW = zeros(m,1);
for i=1:m
    W=cell2mat(W_cell(i));
    u = norm(W(1:d,d+1),2)^2;
    
    
    % 检验 其实放在这里大概率没用, 可以之后再用上.
%     if rank(W) ~= 1 
%        pause;
%     end
    
    [V,D] = eig(W);
    D(D > -0.0001 & D < 0)=0;
    D(D < 0.0001 & D > 0)=0;
    L=V*(D^0.5);
    u1 = trace(L'*L)-1;
    if abs(u-u1)>0.0001
        disp('norm计算出现错误, 两种计算方式不匹配')
        disp(u-u1);
%         pause;
    end
    
    % 此处的u和tr(A'A)-1 一定是相等的.
    norm_AllW(i) = u;
    
    
    
    
end





end