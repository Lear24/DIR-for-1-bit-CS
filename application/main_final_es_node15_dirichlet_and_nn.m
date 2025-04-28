% the all final result for 15-patient empirical study.

% same sample size



%% Dirichlet sample -------------------------------------------------------
%% Channel 1 q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 200;T2 = T1-1+d;Channel_i = 1;
real_beta  = zeros(m,d);

% same sample size
% nn = 300;
% n=nn*ones(m,1);

% Dirichlet
rng(0)
alpha = 0.7;N = 300*m;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-150*m)*r(:,1)/sum(r(:,1)));
n= n + 150;

disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end

rng(1)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_Dirichlet07_Channel1_T200_d40_lam03.mat')

%% Channel 4 q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 440;T2 = T1-1+d;Channel_i = 4;
real_beta  = zeros(m,d);

% same sample size
% nn = 300;
% n=nn*ones(m,1);

% Dirichlet
rng(0)
alpha = 0.7;N = 300*m;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-150*m)*r(:,1)/sum(r(:,1)));
n= n + 150;

disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end

rng(2)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_Dirichlet07_Channel4_T440_d40_lam03.mat')
%% Channel 10 q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 920;T2 = T1-1+d;Channel_i = 10;
real_beta  = zeros(m,d);

% same sample size
% nn = 300;
% n=nn*ones(m,1);

% Dirichlet
rng(0)
alpha = 0.7;N = 300*m;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-150*m)*r(:,1)/sum(r(:,1)));
n= n + 150;

disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end

rng(3)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_Dirichlet07_Channel10_T920_d40_lam03.mat')

%% Channel 1 q=0.95/0.025 with same samples.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 200;T2 = T1-1+d;Channel_i = 1;
real_beta  = zeros(m,d);

% same sample size
nn = 300;
n=nn*ones(m,1);

% Dirichlet
% rng(0)
% alpha = 0.7;N = 300*m;
% r = util_drchrnd([alpha,1-alpha],m);
% n = round((N-150*m)*r(:,1)/sum(r(:,1)));
% n= n + 150;
% 
% disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end

rng(1)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_nn300_Channel1_T200_d40_lam03.mat')
%% Channel 4 q=0.95/0.025 with same sample.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 440;T2 = T1-1+d;Channel_i = 4;
real_beta  = zeros(m,d);

% same sample size
nn = 300;
n=nn*ones(m,1);

% Dirichlet
% rng(0)
% alpha = 0.7;N = 300*m;
% r = util_drchrnd([alpha,1-alpha],m);
% n = round((N-150*m)*r(:,1)/sum(r(:,1)));
% n= n + 150;

% disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end
rng(2)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_nn300_Channel4_T440_d40_lam03.mat')
%% Channel 10 q=0.95/0.025 with same sample.
clear
clc
rng('default')
% generate sample
load('expirical_stady_node15.mat')
m=15;d = 40;T1 = 920;T2 = T1-1+d;Channel_i = 10;
real_beta  = zeros(m,d);

% same sample size
nn = 300;
n=nn*ones(m,1);

% Dirichlet
% rng(0)
% alpha = 0.7;N = 300*m;
% r = util_drchrnd([alpha,1-alpha],m);
% n = round((N-150*m)*r(:,1)/sum(r(:,1)));
% n= n + 150;

disp(num2str(sum(n)/m))


for i = 1:m
    temp = cell2mat(es_data(i));
    real_beta(i,:) = temp(Channel_i,T1:T2);
end
rng(3)
nv=0.3;
% sample generation
[X_cell,Sigma_X,XN] = generation_sample(m,n,d,nv);

% noise generation
ratio2 = 1/2; % the ratio of client with [sigma1,q1]
sigma1 = 0.1; sigma2 = 0.2;
q1 = 0.95;%pisitive
q2 = 0.025;%negtive
[sigma,q] = assign_noise_variance(m,[sigma1,sigma2],[q1,q2],ratio2);% the variance of noise and the probability of sign flipping in each client
[e_cell, xi_cell] = generate_noise(m,n,sigma,q);

% response generation
[Y_cell,YN] = generate_response(m,n,d,real_beta,X_cell,e_cell,xi_cell);

% try all
lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 400; K_max= 10;
[W_eastimated1,beta_estimated1,~]=...
    onebitCS_invex_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda); %(A)

 T_max = 3000; 
[W_eastimated2,beta_estimated2,~]=...
    onebitCS_invexCentering_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,lambda); %(A)Centering
T_max = 400; K_max= 10;
[beta_estimated3,~]=...
    onebitCS_original_Nesterov(beta_initial,X_cell,Y_cell,m,n,d,T_max,K_max,lambda);%(O)
beta_estimated4 = onbitCS_LS(X_cell,Y_cell,m,n,d);
[beta_estimated5,cos_value5,~]=onebitCS_distributed_decoding(X_cell,Y_cell,real_beta,m,n,d,T_max);
[beta_estimated6,cos_value6] = onbitCS_LS_all(XN,YN,real_beta,m);
%         toc

[cos_value1, L2_error1] = util_similarity(beta_estimated1, real_beta,m);
[cos_value2, L2_error2] = util_similarity(beta_estimated2, real_beta,m);
[cos_value3, L2_error3] = util_similarity(beta_estimated3, real_beta,m);
[cos_value4, L2_error4] = util_similarity(beta_estimated4, real_beta,m);
L2_error5 = 2 - 2*abs(cos_value5);
L2_error6 = 2 - 2*abs(cos_value6);
disp(mean(abs([cos_value1,cos_value2,cos_value3,cos_value4,cos_value5,cos_value6])));

cos_value1 = abs(cos_value1);cos_value2 = abs(cos_value2);cos_value3 = abs(cos_value3);
cos_value4 = abs(cos_value4);cos_value5 = abs(cos_value5);cos_value6 = abs(cos_value6);

clear es_data
save('final_test_es_node15_nn300_Channel10_T920_d40_lam03.mat')