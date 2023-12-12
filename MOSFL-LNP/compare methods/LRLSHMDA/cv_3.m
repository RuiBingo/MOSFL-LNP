clc;clear;
tic

%% 原始数据
interaction = importdata('microbedata\interaction.mat');
[nl,nd] = size(interaction);
interaction_ori = interaction;
n = 1;
m=1;
%% 主方法
lambdaM = 1;
lambdaD = 1;
lw = 0.7;
% for lw = 0.1:0.1:0.9
[ KD,KM ] = GIPSim(interaction,1,1 );
sumM=sum(KM,2);
sumD=sum(KD,2);
DM=diag(sumM,0);
DD=diag(sumD,0);
DM_2=diag(sumM.^(-1/2),0);
DD_2=diag(sumD.^(-1/2),0);
LM=DM_2*(DM-KM)*DM_2;
LD=DD_2*(DD-KD)*DD_2;
FM = KM *pinv(KM + lambdaM * (LM*KM))*interaction';
FD = KD *pinv(KD + lambdaD * (LD*KD))*interaction;
F_ori=lw * FM' + (1 - lw) * FD;
%%  5-fold 交叉验证
auc = zeros(1,20);
% rng(0);
for k = 1:20
    k
    lw
    indices = crossvalind('Kfold', nd, 5 );
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(:,index_2) = 0;
       %%%计算得分矩阵
%% 主方法
[ KD,KM ] = GIPSim(interaction,1,1 );
sumM=sum(KM,2);
sumD=sum(KD,2);
DM=diag(sumM,0);
DD=diag(sumD,0);
DM_2=diag(sumM.^(-1/2),0);
DD_2=diag(sumD.^(-1/2),0);
LM=DM_2*(DM-KM)*DM_2;
LD=DD_2*(DD-KD)*DD_2;
FM = KM *pinv(KM + lambdaM * (LM*KM))*interaction';
FD = KD *pinv(KD + lambdaD * (LD*KD))*interaction;
F=lw * FM' + (1 - lw) * FD;
F_ori(:,index_2) = F(:,index_2);
interaction = interaction_ori;
end
%% 画auc曲线 
    pre_label_score = F_ori(:);
    label_y = interaction_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');
%     [ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.95' );
%     figure;
%     aupr(k)=pr_cure(pre_label_score,label_y,'red');
end
%% 
 auc_ave = mean(auc);
 auc_std = std(auc);
 y(m,n)=auc_ave;
 n=n+1;
% end 
 m=m+1;
 n=1;
%  end 
 o=toc;
 o
