function [real_beta, real_theta_para]=generate_random_real_data_rot(m,d,theta_1,theta_2)
% 本函数可以生成生成一个m个d 维度的参数
% 生成方式如下
% 首先生成一个主向量作为主参数 omega_1 (随机生成,每个坐标上都是一个 0-1两点分布)
% 然后生成 m-1 个旋转轴
% 将 参数omega_1关于旋转轴 i = 1,...,m-1 旋转一个角度 tau_i
% tau 是一个随机数,来自(0,theta) 的均匀分布
% 以旋转轴进行schmidt正交化 得到Q_i,找到一个新的空间,将主参数关于这个旋转轴旋转 tau (如果关于关于旋转轴和主参数确定的平面呢?)
% 旋转举证为T=[1,0,0,0...,0;0,cos,-sin,0,...,0;0,0,sin,cos,0,...,0;I
% 参数 omega_i = Q_i'TQ_i omega_1
omega_1=binornd(1,0.5,d,1); % d*1 vector
real_beta=zeros(m,d); % 我们的真实参数矩阵是以 m*d 的矩阵储存的
real_beta(1,:)=omega_1';
I=eye(d);
real_theta_para = zeros(m,1);



for i = 2:m
    error_contol=1;
    rot_axis = binornd(1,0.5,d,1); %旋转轴 - n 维关于旋转轴的旋转不是唯一的 我们可以有不同的确定旋转平面的方式
    tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
    tau = tau*(2*binornd(1,0.5,1,1)-1);
    %     disp(cos(tau));
    %     T=eye(d);
    %     T(2,2)=cos(tau);T(2,3)=-sin(tau);
    %     T(3,2)=sin(tau);T(3,3)=cos(tau);
    %     A=[rot_axis,I];
    A=[omega_1,rot_axis,I];
    %     A=[omega_1,I]; % 这样也可以产生旋转了tau角度的向量,但是总是绕着同一个轴,无法引入更多的随机性. 至少应该随机取d-2个平面作为 旋转轴
    [~,J]=rref(A); % J 中为A的极大线性无关组的坐标信息
    %     disp(rank(A)==d);
    %     disp(J(1)==1);
    Q = util_gram_schmidt(A(:,J)); % 正交化
    real_beta(i,:)= util_rotmnd(Q(:,3:d),tau)*omega_1;
    
    while error_contol==1
        if (real_beta(1,:)*real_beta(i,:)'/norm(real_beta(1,:),2)/norm(real_beta(i,:),2)-cos(tau))<10^(-4)
            fprintf('设备 %d 与原向量角度为 %1.4f Pi \n',i,tau/pi)
            real_theta_para(i) = tau/pi;
            error_contol=0;
        else
            fprintf('设备%d 无法确定旋转角度, 建议重新生成\n',i)
            rot_axis = binornd(1,0.5,d,1);
            tau= theta_1 + (theta_2 - theta_1)*rand(1,1);
            A=[omega_1,rot_axis,I];
            [~,J]=rref(A); % J 中为A的极大线性无关组的坐标信息
            Q = gram_schmidt(A(:,J)); % 正交化
            real_beta(i,:)= util_rotmnd(Q(:,3:d),tau)*omega_1;
        end
        
        
    end 
end








end