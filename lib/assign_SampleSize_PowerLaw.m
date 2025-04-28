function [n,N] = assign_SampleSize_PowerLaw(para_mean,para_var,sample_mean,sample_range_ratio,m)
% n = assign_SampleSize_PowerLaw(para_mean,para_var,sample_mean,sample_range_ratio,m)
% 
% para_mean: the mean of lognorm distribution
% para_var: the var of lognorm distribution
% 
% the parameter of lognormal
% pl_mu = log(para_mean^2/sqrt(para_var+para_mean^2));
% pl_sig = sqrt(log(para_var/(para_mean^2) + 1));
% 
% sample_mean: the expected mean of sample size is  (sample_mean - para_mean);
% sample_range_ratio: Adjust the unbalance of data on clients
% m: number of client.
% 
% 
% 


% Obtain lognormal parameters corresponding to the mean and variance
pl_mu = log(para_mean^2/sqrt(para_var+para_mean^2));
pl_sig = sqrt(log(para_var/(para_mean^2) + 1));
disp([pl_mu,pl_sig])


n = (lognrnd(pl_mu, pl_sig, [m, 1])+sample_mean)*sample_range_ratio;

n = ceil(n); 
N = sum(n); 

sorted_samples = sort(n);

% 绘制排序后的数据
figure(10);
plot(sorted_samples, '-bo', 'LineWidth', 1.5);
grid on;
hold on;
plot([0,m],[floor(N/m),floor(N/m)], '--r', 'LineWidth', 1.5);

legend('Sample Size in Power Law Case', 'N/m', 'Location', 'best')
% 添加图形标题和标签
title('Sorted Samples Per Client', 'FontSize', 14);
xlabel('Clinet Index (Sorted)', 'FontSize', 12);
ylabel('Sample size on each client', 'FontSize', 12);




end 