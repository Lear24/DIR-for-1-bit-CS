function U = generate_transfer_U(W_cell,m,d)
% U = generate_transfer_U(W_cell,m,d)
% In the distributed case, we cannot transfer parameters directly, but we can transfer parameters and values.
% This function is used to generate the parameters and U that we send to the client through the server
% 
% U: the matrix tranfer from server to clinets
% 
% W_cell: (m×1)cell with (d+1 × d+1) matrix
% m: number of client
% d: dimension of coefficient


d1 = d+1;
U = zeros(d1,d1);
[~,I2] = util_matrix_I(d);

for m_i = 1:m
    W = cell2mat(W_cell(m_i));
    tr_i = trace(W) - 1;
    U = U + I2*W*I2/tr_i;

end

end