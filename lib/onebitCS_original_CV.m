function [beta_estimated,Accuracy_rate,lam_information,objective_value]=...
    onebitCS_original_CV(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda)
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
% lambda: given regularization parameter (lam_num × 1)
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
beta_initial1 = beta_initial;
beta_estimated=beta_initial;

% (m×1), record the cost function on each client.
cost_all_client=zeros(m,1);
% (T_max×1), record the cost function on each client.
objective_value = zeros(T_max+1,1);

% assign train set and test set
[X_train,Y_train, X_test,Y_test,n_train,n_test] = assign_train_set(X_cell,Y_cell,m,n,d);

lam_num = max(size(lambda));
lam_ite = 1;
Accuracy_rate = zeros(T_max,1);
lam_information = zeros(lam_num,1);
lam_information(1) = 1;

%The initial objevtive function valiue ------------------------------------
for lam_ite = 1:lam_num
    beta_estimated = beta_initial1;
    lam = lambda(lam_ite);
    for m_i = 1:m
        X_i=cell2mat(X_train(m_i));
        Y_i=cell2mat(Y_train(m_i));
        initial_params =  cell2mat(beta_estimated(m_i));
        cost_all_client(m_i)=objective_original(initial_params,beta_estimated,X_i,Y_i,m,n(m_i),d,m_i,lam);
    end
    objective_value(1) = sum(abs(cost_all_client));
    disp(objective_value(1));
    alpha = 0.0015;
    gamma = 0.9;
    % The conmunication round -------------------------------------------------
    for T=1:T_max
        lam = lambda(lam_ite);
        cost_temp=zeros(m,1);
        beta_initial = beta_estimated;
        % the client part -----------------------------------------------------
        for m_i=1:m
            beta_cell = beta_initial;
            
            X_i=cell2mat(X_train(m_i));
            Y_i=cell2mat(Y_train(m_i));
            n_i = n_train(m_i);
            initial_params = cell2mat(beta_cell(m_i));
            initial_params = initial_params';
            
            % optimalization part on local client------------------------------
            for K = 1:K_max
                
                if K==1
                    [~ ,grad]= objective_original(initial_params',beta_cell,X_i,Y_i,m,n_i,d,m_i,lam);
                    obparams = initial_params - alpha*grad;
                    v = obparams + gamma*(obparams - initial_params);
                else
                    [~ ,grad]= objective_original(v',beta_cell,X_i,Y_i,m,n_i,d,m_i,lam);
                    initial_params = obparams;
                    obparams = v - alpha*grad;
                    v = obparams + gamma*(obparams - initial_params);
                end
                
                
            end
            [J, ~] = objective_original(obparams',beta_cell,X_i,Y_i,m,n_i,d,m_i,lam);
            cost = J;
            % feadback the optimalization result on local client---------------
            beta_estimated(m_i) = mat2cell(obparams',1,d);
            cost_temp(m_i)=cost(end); %record the objective function value on client i
            
        end
        
        beta = cell2mat(beta_estimated);
        Accuracy_rate(T) = util_accuracy(beta, X_test,Y_test,m, n_test);
        
        
        %     if T == 1
        %
        %     elseif abs(Accuracy_rate(T-1) - Accuracy_rate(T) ) < 10^(-6)
        %         if lam_ite == lam_num
        %             fprintf('在轮次%d时, 所有的正则化参数遍历完成\n',T);
        %             lam_information(end) = T;
        %             break;
        %         else
        %             lam_ite = lam_ite + 1;
        %             lam_information(lam_ite)=T;
        %         end
        %     end
        cost_all_client=cost_temp;
        objective_value(T+1) = sum(abs(cost_all_client)); % record objective function value
        
        
        
        
        
        %     % At the server, Broadcast the value of the objective function---------
        %     if  mod(T,100)==0
        %         fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        %     end
        %     % At the server, determine whether to exit the loop--------------------
        %     if abs(sum(cost_temp)-sum(cost_all_client))<10^(-4)
        %         fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        %         fprintf('算法收敛,总共通讯了%d轮\n',T);
        %         % 算法收敛, 记录最终的估计值
        %         objective_value(T+1:end) = sum(abs(cost_temp));% 算法收敛, 记录最终的目标函数值
        %         break;
        %     elseif sum(cost_temp)-sum(cost_all_client)>0.1*sum(cost_all_client)
        %         fprintf('损失函数增加,算法收敛,总共通讯了%d轮\n',T);
        %         fprintf('在轮次%d时, 目标函数值cost为 %f \n',T,sum(cost_temp));
        %         beta_estimated = beta_initial;
        %         objective_value(T+1:end) = sum(abs(cost_temp));
        %         break;
        %     else
        %         % continue loop-----------
        %         cost_all_client=cost_temp;
        %         objective_value(T+1) = sum(abs(cost_all_client)); % record objective function value
        %     end
        %     % The conmunication round ---------------------------------------------
        
    end
end
beta_estimated = beta;
end


function accu_rate =  util_accuracy(beta_mat, X_test, Y_test,m,n_test)

accu_rate = 0;
for i = 1:m
    X = cell2mat(X_test(i));
    Y = cell2mat(Y_test(i));
    beta = beta_mat(i,:);
    nn_test = n_test(i);
    
    
    Y_predict  = sign(X*beta');
    
    accu_rate = accu_rate + norm(Y - Y_predict)/nn_test;
    
end
accu_rate = accu_rate/m;
end

function [X_train,Y_train, X_test,Y_test,n_train,n_test] = assign_train_set(X_cell,Y_cell,m,n,d)

ratio = 0.8;
X_train = cell(m,1);
Y_train = cell(m,1);
X_test  = cell(m,1);
Y_test  = cell(m,1);

n_train =zeros (m,1);
n_test =zeros (m,1);


for i = 1:m
    X = cell2mat(X_cell(i));
    Y = cell2mat(Y_cell(i));
    
    nn = n(i);
    nn_train = round(ratio*nn);
    nn_test = nn - nn_train;
    
    n_train(i) = nn_train;n_test(i) = nn_test;
    
    ind = randsample(nn,nn_train);
    
    X_train(i) = mat2cell( X(ind,:),nn_train,d);
    Y_train(i) = mat2cell( Y(ind,:),nn_train,1);
    
    X(ind,:) = [];Y(ind) = [];
    
    X_test(i) = mat2cell( X,nn_test,d);
    Y_test(i) = mat2cell( Y,nn_test,1);
    
end



end
