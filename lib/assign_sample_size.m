function n = assign_sample_size(N,m,r,ratio)
% n = assign_sample_size(N,m,r,ratio)
% n - m-varibles (m×1)
% N - all sample szie
% m - number of clients
% ratio - the ratio of clients with 2*(1-r)*N/m and 2*r*N/m

n1 = ceil(ratio*m);
n=zeros(m,1);
n(1:(m-2*n1),1) = N/m;
for i = 1:n1
    n((m-2*n1+1):(m-n1),1) = round(2*(1-r)*N/m);
    n((m-n1+1):end,1) = round(2*r*N/m);
end

if sum(n)~= N
   fprintf("样本量分配有误 \n")
   pause;
end

end