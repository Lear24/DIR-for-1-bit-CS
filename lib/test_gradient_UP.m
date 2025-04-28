function  [vdP, dP]= test_gradient_UP(W,U_i,m,d,e)
% [vdP, dP]= test_gradient_UP(W,U_i,m,d,e)



[~,A] = util_W2A(W,d);
dP=zeros(d+1,d+1);


for i=1:d+1
    for j=1:d+1
    I=zeros(d+1,d+1);
    I(i,j)=e;
    A_1 = A+I;
    A_2 = A-I;
    W_1 = A_1*A_1';
    W_2 = A_2*A_2';
    p1=fun_penalty_UP(W_1,U_i,m,d);
    p2=fun_penalty_UP(W_2,U_i,m,d);
    dP(i,j)=(p1-p2)/(2*e);
    end
end

dP(d+1,d+1)=0;
d1 = d+1;
vdP =zeros(d1*d1,1);
for i = 1:d1*d1
    vdP(i) = dP(i);
end



end

