clear
clc

load('final_test_es_patient1_Dirichlet07_eeg3_T1000_d200_lam03.mat')


result_box_cos = zeros(m,6);
result_box_L2 = zeros(m,6);

result_box_cos(:,1) = cos_value1;result_box_cos(:,2) = cos_value2;result_box_cos(:,3)  = cos_value3;
result_box_cos(:,4)  = cos_value4;result_box_cos(:,5) = cos_value5;result_box_cos(:,6)  = cos_value6;

result_box_L2(:,1) = L2_error1;result_box_L2(:,2) = L2_error2;result_box_L2(:,3) = L2_error3;
result_box_L2(:,4) = L2_error4;result_box_L2(:,5) = L2_error5;result_box_L2(:,6) = L2_error6;

result_box_cos(:,3)=[];

figure(1)
% boxplot(result_box_cos, 1:5,'Whisker',1);
boxplot(result_box_cos, 1:5);

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
saveas(1,'test_es_cos_onepatient_eeg3.png')

%L2
result_box_L2(:,3)=[];
figure(2)
% boxplot(result_box_L2, 1:5,'Whisker',1);
boxplot(result_box_L2, 1:5);

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
saveas(2,'test_es_L2_onepatient_eeg3.png')