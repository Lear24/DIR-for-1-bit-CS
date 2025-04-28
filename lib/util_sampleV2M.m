function Q = util_sampleV2M(X,Y,n,d)
% tranfer to sample in form of invex matrix relaxation
% from (d×1) vector to (d+1 × d+1)matrix.
% X: (n*d)
% Y: (n*1)
% already averaged by sample size

Q=zeros(d+1);
for i=1:n
    Q=Q+[X(i,:)' ; -Y(i)]*[X(i,:) , -Y(i)];
end
Q=Q/n;

end