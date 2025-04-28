clear
clc

load('final_test_es_node15_nn300_Channel4_T440_d40_lam03.mat')
result_box_cos = zeros(m,6);
result_box_L2 = zeros(m,6);

result_box_cos(:,1) = cos_value1;result_box_cos(:,2) = cos_value2;result_box_cos(:,3)  = cos_value3;
result_box_cos(:,4)  = cos_value4;result_box_cos(:,5) = cos_value5;result_box_cos(:,6)  = cos_value6;

result_box_L2(:,1) = L2_error1;result_box_L2(:,2) = L2_error2;result_box_L2(:,3) = L2_error3;
result_box_L2(:,4) = L2_error4;result_box_L2(:,5) = L2_error5;result_box_L2(:,6) = L2_error6;

result_box_cos(:,2)=[];
figure(1)

% boxplot(result_box_cos, 1:5,'Whisker',1);
h=boxplot(result_box_cos, 1:5);
% boxplot(result_box_cos, 1:5,'Jitter',0.2);

%如何单独处理异常值

% 获取所有异常值的句柄
outliers = findobj(h, 'Tag', 'Outliers');

% 针对第三个箱线图（即第三列的数据），调整其异常值位置
% outliers 是一个包含所有异常值的句柄数组，其中每个元素对应一个箱线图中的异常值
outliers_data = outliers(3);

% 获取当前异常值的X坐标
outlier_positions = get(outliers_data, 'XData');

% 给异常值添加一定的偏移量 (jitter)
jitter_amount = 0.05;  % 偏移量的大小

rng(1) % 随机错开
a = randn(size(outlier_positions));

outlier_positions = outlier_positions + jitter_amount *a ;
set(outliers_data, 'XData', outlier_positions);

% colors = [0,0.8,1;0,0.8,0;1,1,0.2;0.7,0.1,0.1;0.2,0.4,0.8];
colors = [0,0.8,1;0,0.8,0;1,1,0.2;0.7,0.1,0.1;0.2,0.4,0.8;0,0,1];
h = findobj(gca,'Tag','Box');
f = fliplr(1:5);
for i=1:length(h)
    j = f(i);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

xticks([1 2 3 4 5]);
xticklabels({'DIR', 'CIR', 'SLS', 'DRD', 'PLS'});
% legend('DIR', 'CIR', 'SLS', 'DRD', 'PLS','Location','best')
   yticks(0:0.1:1);
set(gca, 'FontSize', 20);  
saveas(1,'test_es_cos_node15_nn300_channel_4.png')

%L2
result_box_L2(:,3)=[];
figure(2)
% boxplot(result_box_L2, 1:5,'Whisker',1);
h = boxplot(result_box_L2, 1:5);
% boxplot(result_box_L2, 1:5,'Jitter',0.2);

% 获取所有异常值的句柄
outliers = findobj(h, 'Tag', 'Outliers');

% 针对第三个箱线图（即第三列的数据），调整其异常值位置
% outliers 是一个包含所有异常值的句柄数组，其中每个元素对应一个箱线图中的异常值
outliers_data = outliers(3);

% 获取当前异常值的X坐标
outlier_positions = get(outliers_data, 'XData');

% 给异常值添加一定的偏移量 (jitter)
jitter_amount = 0.05;  % 偏移量的大小


outlier_positions = outlier_positions + jitter_amount *a ;
set(outliers_data, 'XData', outlier_positions);


colors = [0,0.8,1;0,0.8,0;1,1,0.2;0.7,0.1,0.1;0.2,0.4,0.8];
% colors = [0,0.8,1;0,0.8,0;1,1,0.2;0.7,0.1,0.1;0.2,0.4,0.8;0,0,1];
h = findobj(gca,'Tag','Box');
f = fliplr(1:5);
for i=1:length(h)
    j = f(i);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

xticks([1 2 3 4 5]);
xticklabels({'DIR', 'CIR', 'SLS', 'DRD', 'PLS'});
% legend('DIR', 'CIR', 'SLS', 'DRD', 'PLS','Location','best')
yticks(0:0.2:4);
set(gca, 'FontSize', 20);
saveas(2,'test_es_L2_node15_nn300_channel_4.png')