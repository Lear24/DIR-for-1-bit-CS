function  dP = test_gradient_vP(beta_all,e,m,d,m_i)
% dP = test_gradient_vP(beta_all,e,m,d,m_i)
% 近似生成梯度
% 使用梯度的定义近似产生梯度
% (f(x+e)-f(x-e))/2e
% 
%
dP=zeros(1,d);
beta=cell2mat(beta_all(m_i));
for i=1:d
    I=zeros(1,d);p1 = 0;p2 = 0;
    beta_all_1=beta_all;
    beta_all_2=beta_all;
    I(1,i)=e;
    beta_all_1(m_i)=mat2cell(beta+I,1);
    beta_all_2(m_i)=mat2cell(beta-I,1);
%     for j = 1:m
%         p1= p1 + fun_penalty_vecP(beta_all_1,m,j);
%         p2= p2 + fun_penalty_vecP(beta_all_2,m,j);
%     end
    p1= p1 + fun_penalty_vecP(beta_all_1,m,m_i);
    p2= p2 + fun_penalty_vecP(beta_all_2,m,m_i);
    dP(1,i)=(p1-p2)/(2*e);
end
dP = dP';
end