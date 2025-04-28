function beta=assign_initialization2d (m,d,theta1,theta2)
% [beta]=initialization(k,d)
% 按照每个位置上生成0-1的方式来进行初始化，这个在离beta距离很远的时候其实不一定适用的
% 这个地方我们采用了和生成real_beta
% beta=mat2cell(binornd(1,0.5,k,d),repelem(1,k));
% theta = (1:m)*pi/(m+1)+rand(1,1)*pi;
theta = (1:m)*pi/(m+1)+rand(1,1)*pi;
% theta = theta1 + (1:m)*(theta2 - theta1)/(m+1)+rand(1,1)*pi;
beta=zeros(m,d);
v=[1;0];
for i  = 1:m
    b = util_rot2d(v, theta(i));
    beta(i,:) = b;
end

end