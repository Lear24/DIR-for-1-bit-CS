function [W,A] = util_A2W(para_A,d)
% [W,L] = util_A2W(para_A,d)
% Obtain the matrix W according to A in the form of vector.
% And record the matrix A as a vector
% 
% W: map to A, W = AA';
% para_A: record the matrix A as a vector (d+1 × d+1)

d1=d+1;
A=zeros(d1,d1);
for i = 1:d1*d1
   A(i) = para_A(i) ;
end
W=A*A';
end