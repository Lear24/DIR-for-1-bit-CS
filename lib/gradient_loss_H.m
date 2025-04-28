function [vdH, dH]= gradient_loss_H(W,X,Y,n,d)
% [vdH, dH]= gradient_loss_H(W,X,Y,n,d)
% The gradient ofloss function in invex problem (A) (average by sample size)
% beta: (d*1( the coefficient in vector form
% X: (n×d)
% Y: (n×1)
% 
% dH: the gradient in form of matrix  (d+1 × d+1)
% vdH: the gradient in form of vector ((d+1)^2 × 1)


d1 = d+1;

Q = util_sampleV2M(X,Y,n,d);

[~,A] = util_W2A(W,d);

dH= 2*Q*A; % 每一步都要注意梯度的计算方式!!!
dH(d1,d1)=0;% 不一定会自然等于0, 反正最后一个元素的梯度必须为0

vdH =zeros(d1*d1,1);
for i = 1:d1*d1
    vdH(i) = dH(i);
end

end