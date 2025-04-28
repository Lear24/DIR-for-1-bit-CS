% the all final result for one-patient empirical study.

%% Dirichlet sample -------------------------------------------------------
%% eeg 1 with q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')

%note 第45 通道有点特殊
load('1_20131027.mat') 
djc_eeg1(45,:)=[];


rng(0)
m=61;d=200;T1=3400;T2=T1+d-1;
N = 1500*m;
alpha = 0.7;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-61*1000)*r(:,1)/sum(r(:,1)));
n = n+1000;
disp(num2str(sum(n)/m))

rng(1)
real_beta = djc_eeg1(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_Dirichlet07_eeg1_T3400_d200_lam03.mat')

%% eeg 3 with q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')

%note 第45 通道有点特殊
load('1_20131027.mat') 
djc_eeg3(6,:)=[];
djc_eeg3(44,:)=[];
djc_eeg3(54,:)=[];

rng(0)
m=59;d=200;T1=1000;T2=T1+d-1;
N = 1500*m;
alpha = 0.7;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-59*1000)*r(:,1)/sum(r(:,1)));
n = n+1000;
disp(num2str(sum(n)/m))

rng(2)
real_beta = djc_eeg3(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_Dirichlet07_eeg3_T1000_d200_lam03.mat')

%% eeg 5 with q=0.95/0.025 with Dirichlet.
clear
clc
rng('default')

%note 第45 通道有点特殊
load('1_20131027.mat') 
djc_eeg5(6,:)=[];
djc_eeg5(44,:)=[];
djc_eeg5(54,:)=[];

rng(0)
m=59;d=200;T1=1200;T2=T1+d-1;
N = 1500*m;
alpha = 0.7;
r = util_drchrnd([alpha,1-alpha],m);
n = round((N-59*1000)*r(:,1)/sum(r(:,1)));
n = n+1000;
disp(num2str(sum(n)/m))

rng(3)
real_beta = djc_eeg5(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_Dirichlet07_eeg5_T1200_d200_lam03.mat')



%% same sample ------------------------------------------------------------
%% eeg 1 with q=0.95/0.025 with same sample
clear
clc

rng('default')
rng(1)

%note 第45 通道有点特殊
load('1_20131027.mat') 
T = 400;
djc_eeg1(45,:)=[];




m=61;n=1500*ones(m,1);d=200;T1=3400;T2=T1+d-1;
real_beta = djc_eeg1(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_nn1500_eeg1_T3400_d200_lam03.mat')

%% eeg 3 with q=0.95/0.025 with same sample
clear
clc

rng('default')
rng(2)

%note 第45 通道有点特殊
load('1_20131027.mat') 
T = 400;
djc_eeg3(6,:)=[];
djc_eeg3(44,:)=[];
djc_eeg3(54,:)=[];

m=59;n=1500*ones(m,1);d=200;T1=1000;T2=T1+d-1;
real_beta = djc_eeg3(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_nn1500_eeg3_T1000_d200_lam03.mat')

%% eeg 5 with q=0.95/0.025 with same sample
clear
clc

rng('default')
rng(3)

%note 第45 通道有点特殊
load('1_20131027.mat') 
T = 400;
djc_eeg5(6,:)=[];
djc_eeg5(44,:)=[];
djc_eeg5(54,:)=[];

m=59;n=1500*ones(m,1);d=200;T1=1000;T2=T1+d-1;
real_beta = djc_eeg5(1:m,T1:T2);
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

%  try all

lambda = 0.3;
beta_initial =assign_initialization(m,d);
T_max = 1000; K_max= 10;
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

%
clear djc*
save('final_test_es_patient1_nn1500_eeg5_T1200_d200_lam03.mat')
