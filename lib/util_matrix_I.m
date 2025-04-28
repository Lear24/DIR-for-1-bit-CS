function  [I1,I2] = util_matrix_I(d)
% [I1,I2] = matrix_I(d)
% Generate two constant matrices that we need to use based on the dimension
% d: generate (d+1)-matrix, I1 and I2, with same defination in paper.
d1 = d+1;
I = zeros(d1,1);
I1 = diag(I);
I1(d1,d1) = 1;
I = ones(d1,1);
I2 = diag(I)-I1;



end