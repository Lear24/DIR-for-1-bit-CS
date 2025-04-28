function util_plot_warmup(Accuracy_rate,lam_information,lambda)

lam_num = max(size(lambda));


figure(1)
hold on 
grid on
box on
% plot(1:lam_information(end),Accuracy_rate(1:lam_information(end)),'LineWidth',2)


% plot(1:lam_information(1), Accuracy_rate(1:lam_information(1)),'LineWidth',2)

for i = 2:lam_num
    hold on 
    plot(lam_information(i-1):lam_information(i),1-Accuracy_rate(lam_information(i-1):lam_information(i)),'LineWidth',2)
end




end