function [W_eastimated,beta_estimated,objective_value]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda)
% [W_eastimated,objective_value,L2_error]=onebitCS_invex_gd(beta_initial,X,Y,m,n,d,T_max,K_max,lambda)
%
% Optimalize (A) in distributed setting.
%
% beta_initial: (m×d) the coefficient in vector form on *all* clients.
% X_cell: (m×1)cell with (n×d)
% Y_cell: (m×1)cell with (n×1)
%
% m: number of clinet
% n: (m×1) vector. The sample size on each clients is different.
% d: dimension.
%
% T_max: The max iteration of conmunication round.
% K_max: The max iteration locally.
%
% lambda: given regularization parameter
%
% W_eastimated: The final estimated coefficient of (W) in form of matrix.
% beta_estimated: The final estimated coefficient of (O) in form of vector.
% objective_value: Record the objective value in each conmunication round.
%
% L2_error: The average L2 error in the conmunication round - in real-parameter version
% Pred_error: The average predict error in the conmunication round - in warming-up version
%


% The initial part --------------------------------------------------------
[~,I2] = util_matrix_I(d);

% Given the initialization of (A)
beta_cell=mat2cell(beta_initial,repelem(1,m));
[W_eastimated,~]=util_O2W_cell(beta_cell,m,d); %record the estimated W
U = generate_transfer_U(W_eastimated,m,d); %record the transmitting vector U

% (m×1), record the cost function on each client.
cost_all_client=zeros(m,1);
% (T_max×1), record the cost function on each client.
objective_value = zeros(T_max+1,1);

%The initial objevtive function valiue ------------------------------------
for i = 1:m
    X_i=cell2mat(X_cell(i));
    Y_i=cell2mat(Y_cell(i));
    
    initial_params = util_W2A(cell2mat(W_eastimated(i)),d);
    [W1,L1] = util_A2W(initial_params,d);
    u = trace(L1'*L1)-1;
    U_i =  U - I2*W1*I2/u;
    
    cost_all_client(i)=objective_invex_A(initial_params,U_i,X_i,Y_i,m,n(i),d,lambda);
end
objective_value(1) = sum(abs(cost_all_client));
disp(objective_value(1));
% The conmunication round -------------------------------------------------
% if lambda > 1.6
%     alpha = 0.0005;
% else
%     alpha = 0.002;
% end
alpha = 0.0001;
% if lambda == 1.5
%     K_max = 1;
% end
gamma = 0.9;
v = zeros ((d+1)^2,m);
initial_params = zeros ((d+1)^2,m);
obparams = zeros ((d+1)^2,m);
for T=1:T_max
    W = W_eastimated;
    cost_temp=zeros(m,1);
    U = generate_transfer_U(W_eastimated,m,d); %U: transmission from server to clients.
    
    % the client part -----------------------------------------------------
    parfor m_i=1:m
        
        
        X_i=cell2mat(X_cell(m_i));
        Y_i=cell2mat(Y_cell(m_i));
        n_i = n(m_i);
        initial_params(:,m_i) = util_W2A(cell2mat(W(m_i)),d);
        [W1,L1] = util_A2W(initial_params(:,m_i),d);
        u = trace(L1'*L1)-1;
        U_i =  U - I2*W1*I2/u;
        
        
        % optimalization part on local client------------------------------
        %         v = zeros((d+1)^2,1);

        for K = 1:K_max
            
            if T==1
                [~ ,grad]= objective_invex_A(initial_params(:,m_i),U_i,X_i,Y_i,m,n_i,d,lambda);
                obparams(:,m_i) = initial_params(:,m_i) - alpha*grad;
                v(:,m_i) = obparams(:,m_i) + gamma*(obparams(:,m_i) - initial_params(:,m_i));
            else
                [~ ,grad]= objective_invex_A(initial_params(:,m_i),U_i,X_i,Y_i,m,n_i,d,lambda);
                initial_params(:,m_i) = obparams(:,m_i);
                obparams(:,m_i) = v(:,m_i) - alpha*grad;
                v(:,m_i) = obparams(:,m_i) + gamma*(obparams(:,m_i) - initial_params(:,m_i));
            end
        end
        [J, ~] = objective_invex_A(obparams(:,m_i),U_i,X_i,Y_i,m,n_i,d,lambda);
        cost = J;
        % feadback the optimalization result on local client---------------
        [W_temp,~] = util_A2W(obparams(:,m_i),d);
        W(m_i)=mat2cell(W_temp,d+1); %record the (W) on client i
        cost_temp(m_i)=cost(end); %record the objective function value on client i
        
    end
    
    
    
    % At the server, Broadcast the value of the objective function---------
    if  mod(T,100)==0
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
    end
    % At the server, determine whether to exit the loop--------------------
    if abs(sum(cost_temp)-sum(cost_all_client))<10^(-5)
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        fprintf('算法收敛,总共通讯了%d轮\n',T);
        W_eastimated = W; % 算法收敛, 记录最终的估计值
        objective_value(T+1:end) = sum(abs(cost_temp));% 算法收敛, 记录最终的目标函数值
        break;
    elseif sum(cost_temp)-sum(cost_all_client)>0.01
        fprintf('损失函数增加,算法收敛,总共通讯了%d轮\n',T);
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        objective_value(T+1:end) = sum(abs(cost_temp));
        break;
    else
        % continue loop-----------
        W_eastimated = W;
        cost_all_client=cost_temp;
        objective_value(T+1) = sum(abs(cost_all_client)); % record objective function value
    end
    % The conmunication round ---------------------------------------------
    
end
beta_estimated = util_Wcell2betamat(W_eastimated,m,d);
end
function beta_estimated = util_Wcell2betamat(W,m,d)

beta_estimated = zeros(m,d);
for  mm = 1:m
    
    WW = cell2mat(W(mm));
    beta_estimated(mm,:) = WW(1:d,d+1);
    
end

end