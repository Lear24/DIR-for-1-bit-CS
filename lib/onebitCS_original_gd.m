function [beta_estimated,objective_value]=...
    onebitCS_original_gd(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda)
% [W_eastimated,objective_value,L2_error]=onebitCS_invex_gd(beta_initial,X,Y,m,n,d,T_max,K_max,lambda)
%
% Optimalize (O) in distributed setting.
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
% beta_estimated: The final estimated coefficient of (O) in form of vector.
% objective_value: Record the objective value in each conmunication round.
%
% L2_error: The average L2 error in the conmunication round - in real-parameter version
% Pred_error: The average predict error in the conmunication round - in warming-up version
%


% The initial part --------------------------------------------------------


% Given the initialization of (A)
beta_initial=mat2cell(beta_initial,repelem(1,m));
beta_estimated=beta_initial;

% (m×1), record the cost function on each client.
cost_all_client=zeros(m,1);
% (T_max×1), record the cost function on each client.
objective_value = zeros(T_max+1,1);

%The initial objevtive function valiue ------------------------------------
for m_i = 1:m
    X_i=cell2mat(X_cell(m_i));
    Y_i=cell2mat(Y_cell(m_i));
    initial_params =  cell2mat(beta_estimated(m_i));
    cost_all_client(m_i)=objective_original(initial_params,beta_estimated,X_i,Y_i,m,n(m_i),d,m_i,lambda);
end
objective_value(1) = sum(abs(cost_all_client));
disp(objective_value(1));
alpha = 0.0025;
gamma = 0.9;
% The conmunication round -------------------------------------------------
for T=1:T_max
    cost_temp=zeros(m,1);
    beta_initial = beta_estimated;
    % the client part -----------------------------------------------------
    for m_i=1:m
        beta_cell = beta_initial;
        
        X_i=cell2mat(X_cell(m_i));
        Y_i=cell2mat(Y_cell(m_i));
        n_i = n(m_i);
        initial_params = cell2mat(beta_cell(m_i));
        initial_params = initial_params';
        
        % optimalization part on local client------------------------------
        for K = 1:K_max
            [~ ,grad]= objective_original(initial_params',beta_cell,X_i,Y_i,m,n_i,d,m_i,lambda);
            obparams = initial_params - alpha*grad;
            initial_params = obparams;
        end
        [J, ~] = objective_original(obparams',beta_cell,X_i,Y_i,m,n_i,d,m_i,lambda);
        cost = J;
        % feadback the optimalization result on local client---------------
        beta_estimated(m_i) = mat2cell(obparams',1,d);
        cost_temp(m_i)=cost(end); %record the objective function value on client i
        
    end
    
    
    
    % At the server, Broadcast the value of the objective function---------
    if  mod(T,100)==0
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
    end
    % At the server, determine whether to exit the loop--------------------
    if abs(sum(cost_temp)-sum(cost_all_client))<10^(-4)
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        fprintf('算法收敛,总共通讯了%d轮\n',T);
        % 算法收敛, 记录最终的估计值
        objective_value(T+1:end) = sum(abs(cost_temp));% 算法收敛, 记录最终的目标函数值
        break;
    elseif sum(cost_temp)-sum(cost_all_client)>0.1*sum(cost_all_client)
        fprintf('损失函数增加,算法收敛,总共通讯了%d轮\n',T);
        fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        beta_estimated = beta_initial;
        objective_value(T+1:end) = sum(abs(cost_temp));
        break;
    else
        % continue loop-----------
        cost_all_client=cost_temp;
        objective_value(T+1) = sum(abs(cost_all_client)); % record objective function value
    end
    % The conmunication round ---------------------------------------------
    

     
    
end
    beta_estimated = cell2mat(beta_estimated);
end
