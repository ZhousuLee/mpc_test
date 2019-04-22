clc
clear all

%% 采样时间、预测时域、控制时域
sample_time=0.05;
N=500;
Np=10;
Nc=10;
Row=5;
%% 汽车参数
L=2.6
%% 参考轨迹，圆形轨迹，分别为x，y，theta
ref_time=zeros(N,1);
ref_traj=zeros(N,3);
ref_vel=zeros(N,2);
y_ref=zeros(3*Np,1);

for i=1:1:N
    ref_time(i,1)=(i-1)*sample_time;
    ref_traj(i,1)=25*sin(0.2*ref_time(i,1));
    ref_traj(i,2)=25+10-25*cos(0.2*ref_time(i,1));
    ref_traj(i,3)=0.2*ref_time(i,1);
    ref_vel(i,1)=5;
    ref_vel(i,2)=0.103627;

end
figure(1)
hold on
plot( ref_traj(:,1), ref_traj(:,2))


%% 真实轨迹
real_traj=zeros(N,3);
real_control=zeros(N,2);
real_traj(1,:)=[-1 5 0];
%% 预测轨迹
predict_traj=zeros(N,3*Np);
predict_control=zeros(N+1,2);
predict_delta_contrl=zeros(N+1,2);

%% 轨迹误差
bias_traj_real_ref=zeros(N,3);



 Q=100*eye(3*Np,3*Np);    
 R=5*eye(2*Nc);
for i=1:1:N
    vel=ref_vel(i,1);
    delta=ref_vel(i,2);
    phi=ref_traj(i,3);
    T=sample_time;
    A_k=[1    0   -vel*sin(phi)*T;
       0    1   vel*cos(phi)*T;
       0    0   1;];
    B_k=[cos(phi)*T   0;
       sin(phi)*T   0;
       tan(delta)*T/L      vel*T/(cos(delta)^2);];
    A_cell=cell(2,2);
    B_cell=cell(2,1);
    A_cell{1,1}=A_k;
    A_cell{1,2}=B_k;
    A_cell{2,1}=zeros(2,3);
    A_cell{2,2}=eye(2);
    B_cell{1,1}=B_k;
    B_cell{2,1}=eye(2);
    A=cell2mat(A_cell);
    B=cell2mat(B_cell);
   C=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;];

    PHI_cell=cell(Np,1);
    THETA_cell=cell(Np,Nc);   
    for j=1:1:Np
        PHI_cell{j,1}=C*A^j;
        for k=1:1:Nc
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B;
            else 
                THETA_cell{j,k}=zeros(3,2);
            end
        end
    end
    PHI=cell2mat(PHI_cell);%size(PHI)=[Nx*Np Nx+Nu]
    THETA=cell2mat(THETA_cell);%size(THETA)=[Nx*Np Nu*(Nc+1)]
    H_cell=cell(2,2);
    H_cell{1,1}=THETA'*Q*THETA+R;
    H_cell{1,2}=zeros(2*Nc,1);
    H_cell{2,1}=zeros(1,2*Nc);
    H_cell{2,2}=Row;
    H=cell2mat(H_cell);
    kesi=zeros(5,1);
    kesi(1:3,1)=real_traj(i,:)'-ref_traj(i,:)';
    if i==1
        kesi(4:5,1)=[0;0];
    else
         kesi(4:5,1)=real_control((i-1),:)'-ref_vel((i-1),:)';
    end
    for k=0:1:Np-1
        if(i+k>N)
            y_ref(3*k+1,1)=ref_traj(N,1)-ref_traj(i,1);
            y_ref(3*k+2,1)=ref_traj(N,2)-ref_traj(i,2);
            y_ref(3*k+3,1)=ref_traj(N,3)-ref_traj(i,3);
        else
            y_ref(3*k+1,1)=ref_traj(i+k,1)-ref_traj(i,1);
            y_ref(3*k+2,1)=ref_traj(i+k,2)-ref_traj(i,2);
            y_ref(3*k+3,1)=ref_traj(i+k,3)-ref_traj(i,3);
        end
    end
    
     error=PHI*kesi-y_ref;
    f_cell=cell(1,2);
    f_cell{1,1}=2*error'*Q*THETA;
    f_cell{1,2}=0;
    f=cell2mat(f_cell);
    
 %% 约束条件
 
 A_t=tril(ones(Nc));
 A_I=kron(A_t,eye(2));
 if i==1
     Ut=kron(ones(Nc,1),[0;0]);
 else
      Ut=kron(ones(Nc,1),real_control(i-1,:)');
 end

    umin=[-15;-1;];%维数与控制变量的个数相同
    umax=[15;1];
    delta_umin=[-0.2;-0.2;];
    delta_umax=[0.2;0.2];
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);
    A_cons_cell={A_I zeros(2*Nc,1);-A_I zeros(2*Nc,1)};

    b_cons_cell={Umax-Ut;-Umin+Ut};
    A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围
    b_cons=cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值 
    
       % 状态量约束
    M=10; 
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);
    lb=[delta_Umin;0];%（求解方程）状态量下界，包含控制时域内控制增量和松弛因子
    ub=[delta_Umax;M];%（求解方程）状态量上界，包含控制时域内控制增量和松弛因子
    
    %% 开始求解过程
      options = optimset('Algorithm','active-set');
%     options = optimset('Algorithm','interior-point-convex'); 
    [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
    
    predict_delta_contrl(i,1)=X(1);
    predict_delta_contrl(i,2)=X(2);
    if i==1
        predict_control(i,:)=0+predict_delta_contrl(i,:);
    else
        predict_control(i,:)=predict_control(i-1,:)+predict_delta_contrl(i,:);
    end
    
    real_control(i,:)= predict_control(i,:)+ref_vel(i,:);
    

    vd11=real_control(i,1);
    vd22=tan(real_control(i,2))/L*vd11;
%     vd11=ref_vel(1,1);vd22=tan(ref_vel(1,2))/L*vd11;
    X00=real_traj(i,:);
    XOUT=dsolve('Dx-vd11*cos(z)=0','Dy-vd11*sin(z)=0','Dz-vd22=0','x(0)=X00(1)','y(0)=X00(2)','z(0)=X00(3)');
    t=sample_time;
    real_traj(i+1,1)=eval(XOUT.x);
    real_traj(i+1,2)=eval(XOUT.y);
    real_traj(i+1,3)=eval(XOUT.z);
    
    
end

plot(real_traj(:,1),real_traj(:,2),'o');
% plot(real_traj(1:6,1),real_traj(1:6,2),'o');

figure(2)
plot(ref_time(:,1),real_control(:,1))

figure(3)
plot(ref_time(:,1),real_control(:,2))













