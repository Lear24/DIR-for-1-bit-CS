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
%%
rng('default')
rng(0)
rep_num = 100; % 实验重复次数

d=20; % beta 维数
n_vary = [50;60;70;80;90;100;110;120]; % 8
% n_vary = [60;80;100;120;140;160;180;200]; % 8
n_num = 8;

m = 30;
theta_1=0;
theta_2=pi/8;
[real_beta,real_theta_para ]= generate_random_real_data_rot(m,d,theta_1,theta_2); %随机生成正向或者反向的
real_beta_cell = mat2cell(real_beta,repelem(1,m));

nv=0.3;
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.75;%pisitive
q2 = 0.125;%negtive

beta_estimated_cell = cell (n_num,rep_num,6);
result = cell (n_num,rep_num,6,2);

lam = [0.6;0.8;1;1.2;1.6;1.8;2.2];
%%

for rep = 11:rep_num
    for i = 1:n_num
        n = n_vary(i)*ones(m,1);
        fprintf('第%d次重复实验, setting of sample size(n) = %d\n',rep,n_vary(i));
        rng(rep)
        [sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
        [X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);
        [e_cell, xi_cell] = generate_noise(m,n,sigma,q);
        [Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);
        
        T_max = 4000; K_max= 10;
        beta_initial =assign_initialization(m,d);
        
        if i == 1
            lambda = 0.8;
        elseif i == 2
            lambda = 1;
        elseif i == 3
            lambda = 1.2;
        elseif i == 4
            lambda = 1.5;    
        elseif i == 5
            lambda = 1.8;    
        elseif i == 6
            lambda = 2.1;
        elseif i == 7
            lambda = 2.4;
        elseif i == 8
            lambda = 3;
        end
        
        
        % optimalization part ---------------------------------------------
        [W_eastimated1,beta_estimated1,objective_value1]=...
            onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)
        [W_eastimated2,beta_estimated2,objective_value2]=...
            onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
        [beta_estimated3,objective_value3]=...
            onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
        beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
        [beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
        [beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
        
        beta_estimated_cell(i,rep,1)=mat2cell(beta_estimated1,m,d);
        beta_estimated_cell(i,rep,2)=mat2cell(beta_estimated2,m,d);
        beta_estimated_cell(i,rep,3)=mat2cell(beta_estimated3,m,d);
        beta_estimated_cell(i,rep,4)=mat2cell(beta_estimated4,m,d);
        beta_estimated_cell(i,rep,5)=mat2cell(beta_estimated5',1,d);
        beta_estimated_cell(i,rep,6)=mat2cell(beta_estimated6',1,d);
        
        [cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
        [cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
        [cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
        [cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
        L2_error5 = 2 - 2*abs(cos_value5);
        L2_error6 = 2 - 2*abs(cos_value6);
        disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));
        
        result(i,rep,1,1) = mat2cell(cos_value1,m,1);result(i,rep,1,2) = mat2cell(L2_error1,m,1);
        result(i,rep,2,1) = mat2cell(cos_value2,m,1);result(i,rep,2,2) = mat2cell(L2_error2,m,1);
        result(i,rep,3,1) = mat2cell(cos_value3,m,1);result(i,rep,3,2) = mat2cell(L2_error3,m,1);
        result(i,rep,4,1) = mat2cell(cos_value4,m,1);result(i,rep,4,2) = mat2cell(L2_error4,m,1);
        result(i,rep,5,1) = mat2cell(cos_value5,m,1);result(i,rep,5,2) = mat2cell(L2_error5,m,1);
        result(i,rep,6,1) = mat2cell(cos_value6,m,1);result(i,rep,6,2) = mat2cell(L2_error6,m,1);
        
        save('test_total_sample_vary_m30_pi8.mat')
    end
end

%%
load('test_total_sample_vary_m30_pi8.mat')
rep_num = 100;
beta_quantail = zeros(n_num,2,4,rep_num);
for i =1 :n_num
    for rep = 1:rep_num
        for k = 1:3 % 1/4 1/2 3/4 quantail;
            
            beta_quantail(i,1,k,rep) = quantile(abs(cell2mat(result(i,rep,1,1))),k/4);
            beta_quantail(i,2,k,rep) = quantile(abs(cell2mat(result(i,rep,4,1))),k/4);
            
        end
        
        beta_quantail(i,2,4,rep) = sum(abs(cell2mat(result(i,rep,1,1))) > abs(cell2mat(result(i,rep,4,1))))/m; %百分比
    end
end
beta_quantail_mean = mean(beta_quantail,4);
beta_quantail_var = var(beta_quantail,1,4);


result_mat_cos = zeros(n_num,6);
result_mat_L2 = zeros(n_num,6);
for i = 1:n_num
    
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
result_mat_L2 = sqrt(result_mat_L2);

figure(1)
hold on
grid on
box on


xlim([min((n_vary)),max((n_vary))])
xticks(n_vary)
ylim([min(result_mat_cos(:,5))-0.1,1])
% plot((n_vary),result_mat_cos(:,1),'-ob','LineWidth',3)
% plot((n_vary),result_mat_cos(:,2),'--*r','LineWidth',1.5)
% % plot((n_vary),result_mat_cos(:,3),'--+g','LineWidth',1)
% plot((n_vary),result_mat_cos(:,4),'-oc','LineWidth',2,'Color','#4DBEEE')
% plot((n_vary),result_mat_cos(:,5),'--g','LineWidth',2)
% plot((n_vary),result_mat_cos(:,6),'-oc','LineWidth',2)

plot(n_vary, result_mat_cos(:,1), 'b-*','LineWidth',2.5, 'MarkerSize', 10);  % 实心星形（红色实线）
hold on;
plot(n_vary, result_mat_cos(:,2), 'r-.d', 'LineWidth',1.5,'MarkerSize', 8);  % 蓝色虚线圆圈
plot(n_vary, result_mat_cos(:,4), 'g--o', 'LineWidth',1.5);  % 绿色点划线三角形
plot(n_vary, result_mat_cos(:,5), 'k-.s', 'Marker', 's', 'MarkerFaceColor', 'k');  % 黑色点划线方形
plot(n_vary, result_mat_cos(:,6), 'm--v', 'Marker', 'v', 'MarkerFaceColor', 'm');  % 品红色点划线菱形

% title('The total sample size', 'FontSize', 14);
% xlabel('The number of samples on each node', 'FontSize', 12);
% ylabel('L2 error', 'FontSize', 12);


% legend('Distributed Invex method','Concentrated Invex method',...
%     'Least square method',...
%     'Distributed Decoding','LS for one','Location','best','FontSize', 9,'Position',[0.7,0.3,0.1,0.1])
% legend('Distributed Invex method','Concentrated Invex method',...
%     'Distributed original method','Least square method',...
%     'Distributed Decoding','LS for one','Location','best')
yticks(0:0.1:1);
set(gca, 'FontSize', 20); 
legend('DIR','CIR','SLS','DRD','PLS',...
    'Location','best','Orientation','horizontal','FontSize', 12,'NumColumns', 3)
saveas(1,'test_cos_total_sample_vary_m30_pi8.png')
% saveas(1,'test_cos_total_sample_vary_m30.pdf')




figure(2)
hold on
grid on
box on
xlim([min((n_vary)),max((n_vary))])
xticks(n_vary)
ylim([0,max(result_mat_L2(:,5))+0.1])

% plot((n_vary),result_mat_L2(:,1),'-ob','LineWidth',3)
% plot((n_vary),result_mat_L2(:,2),'-*r','LineWidth',1.5)
% % plot((n_vary),result_mat_L2(:,3),'-+g','LineWidth',1)
% plot((n_vary),result_mat_L2(:,4),'-oc','LineWidth',2,'Color','#4DBEEE')
% plot((n_vary),result_mat_L2(:,5),'-oc','LineWidth',2)
% plot((n_vary),result_mat_L2(:,6),'-oc','LineWidth',2)

plot(n_vary, result_mat_L2(:,1), 'b-*','LineWidth',2.5,'Marker', '*','MarkerSize', 10);  % 实心星形（红色实线）
hold on;
plot(n_vary, result_mat_L2(:,2), 'r-.d', 'LineWidth',1.5,'MarkerSize', 8);  % 蓝色虚线圆圈
plot(n_vary, result_mat_L2(:,4), 'g--o', 'LineWidth',1.5);  % 绿色点划线三角形
plot(n_vary, result_mat_L2(:,5), 'k-.s', 'Marker', 's', 'MarkerFaceColor', 'k', 'DisplayName', 'DRD');  % 黑色点划线方形
plot(n_vary, result_mat_L2(:,6), 'm--v', 'Marker', 'v', 'MarkerFaceColor', 'm', 'DisplayName', 'CIR');  % 品红色点划线菱形


% title('The Heterogeneity of local sample size', 'FontSize', 14);
% xlabel('The Setting of client number', 'FontSize', 12);
% ylabel('L2 error', 'FontSize', 12);


% legend('Distributed Invex method','Concentrated Invex method',...
%     'Distributed original method','Least square method',...
%     'Distributed Decoding','LS for one','Location','best')
yticks(0:0.2:4);
set(gca, 'FontSize', 20);
legend('DIR','CIR','SLS','DRD','PLS',...
    'Location','best','Orientation','horizontal','FontSize', 12,'NumColumns', 3)
saveas(2,'test_L2_total_sample_vary_m30_pi8.png')
% saveas(2,'test_L2_total_sample_vary_m30_pi8.pdf')