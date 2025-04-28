% read me
%
% 变量设置-最简单的情形
% m - number of clinets
% n - number of samples in one clinet
% d - the dimension of coefficient
%
% real_beta - real coefficients
% estimated_beta - estimated coefficients
% X - samples of observation
% Y - response variables
% X_cell/Y_cell - X and Y  in the form of Cell
%
% xi - the samples of sign-flip variables
% q  - the probability of sign flipping
% e - the samples of noise
% sigma - the variance of noise
%
% 本文档的所有向量都以**行向量**的形式储存在矩阵当中.
%%
clear;clc
addpath(genpath('.\lib'));
%%
rng('default')
rng(0)
rep_num = 100; % 实验重复次数

d=20; % beta 维数
m=30;


theta_1=0;
theta_2=pi/8;
[real_beta,real_theta_para ]= generate_random_real_data_rot(m,d,theta_1,theta_2); %随机生成正向或者反向的
real_beta_cell = mat2cell(real_beta,repelem(1,m));

nv=0.3;
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1;
q1 = 0.75;%pisitive
q2 = 0.125;%negtive
sigma2 = 0.2;
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client

lam = [0.3;0.6;0.9;1.2;1.5;1.8];
lam_num = max(size(lam));


beta_estimated_cell = cell (lam_num,rep_num,4);
result = cell (lam_num,rep_num,4,2);
%%

rng(2)
% N = 2400;r=0.4;ratio1=1/10;
% n = assign_sample_size(N,m,r,ratio1);

% para_mean = 40; para_var = 800;sample_mean=60; sample_range_ratio = 0.82;
% [n,N] = assign_SampleSize_PowerLaw(para_mean,para_var,sample_mean,sample_range_ratio,m);

pl_mu = 4;pl_sig=1.8;sample_mean=40; sample_range_ratio = 0.805; % more unbalance in vasualization
% pl_mu = 4;pl_sig=1;sample_mean=60; sample_range_ratio = 0.74;
[n,N] = assign_SampleSize_PowerLaw_pl(pl_mu,pl_sig,sample_mean,sample_range_ratio,m);
disp(N);
n(end) = n(end) - 1;
N = sum(n);
disp(N);

% saveas(10,'test_data_similarity_lam_rep.pdf')
% saveas(10,'test_data_similarity_lam_rep.png')

%%
for rep = 72:100
    
    rng(rep)

    [X_cell,Sigma_X] = generation_sample(m,n,d,nv);
    [e_cell, xi_cell] = generate_noise(m,n,sigma,q);
    Y_cell = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell); 
    beta_initial =assign_initialization(m,d);
    
    for i = 1:lam_num       
        lambda = lam(i);
        fprintf('第%d次重复实验, setting lambda = %d\n',rep,lambda);
        
        T_max = 4000; K_max= 10;


        
        % optimalization part ---------------------------------------------
        [W_eastimated1,beta_estimated1,objective_value1]=...
            onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)
        [W_eastimated2,beta_estimated2,objective_value2]=...
            onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A) Centering
        [beta_estimated3,objective_value3]=...
            onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
        beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
        
        beta_estimated_cell(i,rep,1)=mat2cell(beta_estimated1,m,d);
        beta_estimated_cell(i,rep,2)=mat2cell(beta_estimated2,m,d);
        beta_estimated_cell(i,rep,3)=mat2cell(beta_estimated3,m,d);
        beta_estimated_cell(i,rep,4)=mat2cell(beta_estimated4,m,d);
        
        [cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
        [cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
        [cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
        [cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
        disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4])));
%         disp(mean(abs([cos_value1,cos_value2,cos_value4])));
        
        result(i,rep,1,1) = mat2cell(cos_value1,m,1);result(i,rep,1,2) = mat2cell(L2_error1,m,1);
        result(i,rep,2,1) = mat2cell(cos_value2,m,1);result(i,rep,2,2) = mat2cell(L2_error2,m,1);
        result(i,rep,3,1) = mat2cell(cos_value3,m,1);result(i,rep,3,2) = mat2cell(L2_error3,m,1);
        result(i,rep,4,1) = mat2cell(cos_value4,m,1);result(i,rep,4,2) = mat2cell(L2_error4,m,1);
        
        save('test_similarity_lambda_vary_pi8_rep100.mat')
    end
end

%%

load('test_similarity_lambda_vary_pi8_rep100.mat')

beta_estimated_cell1 = cell (lam_num,rep_num,6);
result1 = cell (lam_num,rep_num,6,2);

for i =1:4
result1(:,:,i,:) = result(:,:,i,:);
end

for i = 1:lam_num
    
    lambda = lam(i);

    for rep = 1:100
        rng(rep)
        
        [X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);
        [e_cell, xi_cell] = generate_noise(m,n,sigma,q);
        [Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);
        
        T_max = 4000; K_max= 10;
        
        
        beta_initial =assign_initialization(m,d);

        
        % optimalization part ---------------------------------------------
        [beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
        [beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
        
        beta_estimated_cell1(i,rep,5)=mat2cell(beta_estimated5',1,d);
        beta_estimated_cell1(i,rep,6)=mat2cell(beta_estimated6',1,d);
        
        L2_error5 = 2 - 2*abs(cos_value5);
        L2_error6 = 2 - 2*abs(cos_value6);
        disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));
        
        result1(i,rep,5,1) = mat2cell(cos_value5,m,1);result1(i,rep,5,2) = mat2cell(L2_error5,m,1);
        result1(i,rep,6,1) = mat2cell(cos_value6,m,1);result1(i,rep,6,2) = mat2cell(L2_error6,m,1);
%         save('test_similarity_lambda_vary_pi3_rep100_complete.mat')
    end
end


%%
result = result1;
rep_num = 100;
beta_quantail = zeros(lam_num,2,4,rep_num);
for i =1 :lam_num
    for rep = 1:rep_num
        for k = 1:3 % 1/4 1/2 3/4 quantail;
            
            beta_quantail(i,1,k,rep) = quantile(abs(cell2mat(result(i,rep,1,1))),k/4);
            beta_quantail(i,2,k,rep) = quantile(abs(cell2mat(result(i,rep,4,1))),k/4);
            
        end
        
        beta_quantail(i,2,4,rep) = sum(abs(cell2mat(result(i,rep,1,1))) > abs(cell2mat(result(i,rep,4,1))))/m;
    end
end
beta_quantail_mean = mean(beta_quantail,4);
beta_quantail_var = var(beta_quantail,1,4);

%%
result_mat_cos = zeros(lam_num,6);
result_mat_L2 = zeros(lam_num,6);
for i = 1:lam_num
    
    for rep = 1:rep_num
        
        result_mat_cos(i,1) = result_mat_cos(i,1) + mean(abs(cell2mat(result(i,rep,1,1))));
        result_mat_L2(i,1) = result_mat_L2(i,1) + mean(abs(cell2mat(result(i,rep,1,2))));
        result_mat_cos(i,2) = result_mat_cos(i,2) + mean(abs(cell2mat(result(i,rep,2,1))));
        result_mat_L2(i,2) = result_mat_L2(i,2) + mean(abs(cell2mat(result(i,rep,2,2))));
        result_mat_cos(i,3) = result_mat_cos(i,3) + mean(abs(cell2mat(result(i,rep,3,1))));
        result_mat_L2(i,3) = result_mat_L2(i,3) + mean(abs(cell2mat(result(i,rep,3,2))));
        result_mat_cos(i,4) = result_mat_cos(i,4) + mean(abs(cell2mat(result(i,rep,4,1))));
        result_mat_L2(i,4) = result_mat_L2(i,4)  + mean(abs(cell2mat(result(i,rep,4,2))));
        result_mat_cos(i,5) = result_mat_cos(i,5) + mean(abs(cell2mat(result(i,rep,5,1))));
        result_mat_L2(i,5) = result_mat_L2(i,5)  + mean(abs(cell2mat(result(i,rep,5,2))));
        result_mat_cos(i,6) = result_mat_cos(i,6) + mean(abs(cell2mat(result(i,rep,6,1))));
        result_mat_L2(i,6) = result_mat_L2(i,6)  + mean(abs(cell2mat(result(i,rep,6,2))));
        
    end
end

result_mat_cos = result_mat_cos/rep_num;
result_mat_L2 = result_mat_L2/rep_num;
%%
% load('test_similarity_lambda_vary_pi8_compelete_rep100.mat')
figure(1)
hold on
grid on
box on
xlim([min(lam),max(lam)])
xticks(lam)
ylim([min(result_mat_cos(:,5))-0.1,1])
plot(lam,result_mat_cos(:,1),'-ob','LineWidth',3)
plot(lam,result_mat_cos(:,2),'--*r','LineWidth',1.5)
% plot(lam,result_mat_cos(:,3),'--+g','LineWidth',1)
plot(lam,result_mat_cos(:,4),'-o','LineWidth',2,'Color','#4DBEEE')
plot(lam,result_mat_cos(:,5),'--g','LineWidth',2)
plot(lam,result_mat_cos(:,6),'-oc','LineWidth',2)



title('The Heterogeneity of Parameter (\theta = \pi/8)', 'FontSize', 14);
xlabel('Regularization parameter Settings(\lambda)', 'FontSize', 12);
ylabel('Cosine value', 'FontSize', 12);


legend('Distributed Invex method','Concentrated Invex method',...
    'Least square method',...
    'Distributed Decoding','LS for one','Location','best')

% legend('Distributed Invex method','Concentrated Invex method',...
%     'Distributed original method','Least square method',...
%     'Distributed Decoding','LS for one','Location','best')


%%
% saveas(1,'test_cos_similarity_lam_pi8_rep100.png')
% save('test_similarity_lambda_vary_pi8_compelete_rep100.mat')

%%  箱线图 - 数据处理
result_box  = zeros(m,6,lam_num);

for i = 1:lam_num
    
    for rep = 1:rep_num 
        result_box(:,1,i) = result_box(:,1,i) + (abs(cell2mat(result(i,rep,1,1))));  
        result_box(:,2,i) = result_box(:,2,i) + (abs(cell2mat(result(i,rep,2,1))));  
        result_box(:,3,i) = result_box(:,3,i) + (abs(cell2mat(result(i,rep,3,1))));  
        result_box(:,4,i) = result_box(:,4,i) + (abs(cell2mat(result(i,rep,4,1))));  
        result_box(:,5,i) = result_box(:,5,i) + (abs(cell2mat(result(i,rep,5,1))));  
        result_box(:,6,i) = result_box(:,6,i) + (abs(cell2mat(result(i,rep,6,1))));  

    end
    
end
result_box = result_box/rep;
%% 箱线图 - 绘制
result_box_cos = result_box;
result_box_cos(:,3,:)=[];
figure(2)
hold on
grid on
box on
i = 5;
% rng(10)
boxplot(result_box_cos(:,:,i), 1:5);



colors = [0,0.8,1;0,0.8,0;1,1,0.2;0.7,0.1,0.1;0.2,0.4,0.8];
h = findobj(gca,'Tag','Box');
f = fliplr(1:5);
for i=1:length(h)
    j = f(i);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

legend('Distributed Invex method','Concentrated Invex method',...
    'Least square method',...
    'Distributed Decoding','LS for one','Location','best')

%% 
saveas(2,'test_boxchart_similarity_lam_pi8_rep100.png')

