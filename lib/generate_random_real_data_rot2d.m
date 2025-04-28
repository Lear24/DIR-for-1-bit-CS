function [real_beta, real_theta_para]=generate_random_real_data_rot2d(m,d,theta_1,theta_2)

% first: random generation for first client

omega_1=[1;0]; % d*1 vector
real_beta=zeros(m,d); % 我们的真实参数矩阵是以 m*d 的矩阵储存的
real_beta(1,:)=omega_1';
I=eye(d);
real_theta_para = zeros(m,1);
for i = 2:m
    theta = theta_1 + (theta_2 - theta_1)*rand(1,1);
    theta = theta * (2*binornd(1,0.5,1,1)-1);
    if binornd(1,0.5,1,1) == 1
        if theta > 0
            theta = pi - theta;
        else
            theta = -pi - theta;
        end
    end
    
    omega = util_rot2d(omega_1,theta);
    
    
    real_theta_para(i) = theta;
    real_beta(i,:)=omega';
    fprintf('设备 %d 与原向量角度为 %1.4f Pi \n',i,theta/pi)
    
end



end