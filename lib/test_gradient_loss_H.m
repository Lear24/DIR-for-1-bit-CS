function  [vdH, dH]= test_gradient_loss_H(W,X,Y,n,d,e)
% [vdH, dH]= test_gradient_loss_H(W,X,Y,n,d,e)



[~,A] = util_W2A(W,d);
dH=zeros(d+1,d+1);


for i=1:d+1
    for j=1:d+1
    I=zeros(d+1,d+1);
    I(i,j)=e;
    A_1 = A+I;
    A_2 = A-I;
    W_1 = A_1*A_1';
    W_2 = A_2*A_2';
    p1=fun_loss_H(W_1,X,Y,n,d);
    p2=fun_loss_H(W_2,X,Y,n,d);
    dH(i,j)=(p1-p2)/(2*e);
    end
end
d1 = d+1;
dH(d+1,d+1)=0;

vdH =zeros(d1*d1,1);
for i = 1:d1*d1
    vdH(i) = dH(i);
end



end

