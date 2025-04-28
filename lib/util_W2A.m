function [para_A,A] = util_W2A(W,d)
% [para_A,L] = util_W2A(W,d)
% Obtain the matrix A according to W.
% And record the matrix A as a vector
% 
% W: map to A, W = AA';
% para_A: record the matrix A as a vector (d+1 × d+1)

% 本函数旨在将正定矩阵分解成两个下三角矩阵的乘积
% para_vector 是L的向量化之后的形式
d1=d+1;
[V,D] = eig(W);
% 我们把接近于0的部分严格的设置为0, 防止还有符号的影响,主要是避免出现虚数
D(D > -0.0001 & D < 0)=0;
D(D < 0.0001 & D > 0)=0;

% 因为不再存在下三角的问题了 所以这个部分很简洁
A = V*(D^0.5);
para_A = zeros(d1*d1,1);
for i = 1:d1*d1
    para_A(i) = A(i);
end

end