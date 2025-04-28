function ome=generate_transfer_vector(beta_all,m,d)
% ome=generate_transfer_vector(beta_all,d,m)
% In the original problem, with distributed setting
% generate the transmitting vector form server to clients.
%
% beta_all: (m×1) cell with (1×d) vector

ome=zeros(1,d);
for i=1:m
    beta = cell2mat(beta_all(i));
    ome=ome+ beta'*beta/norm(beta)^2;
end
end