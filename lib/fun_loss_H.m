function H = fun_loss_H(W,X,Y,n,d)
% F = fun_loss_H(W,X,Y,n,d)
% The loss function in invex problem (W) (average by 2 times of client number)
% To obtain the result of (A), transfer A to W = AA^T.
% W: the coeffcient in form of invex matrix relaxation. (d+1 × d+1)
% X: (n*d)
% Y: (n*1)


Q = util_sampleV2M(X,Y,n,d);

H = trace(Q*W);


end