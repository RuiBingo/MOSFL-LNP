clc;clear;
tic

%% 原始数据
interaction = importdata('microbedata\interaction.mat');
[dn,dr] = size(interaction);
interaction_ori = interaction;
n = 1;
m=1;
beta = 0.01;
%% 主方法
[ KD,KM ] = GIPSim(interaction,1,1 );
Sk2 = beta * interaction' + beta^2 * (KM*interaction' + interaction'* KD);
F_ori=Sk2';
%%  5-fold 交叉验证
auc = zeros(1,20);
% rng(0);
for k = 1:20
    k
    indices = crossvalind('Kfold', dn, 5 );
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(index_2,:) = 0;
       %%%计算得分矩阵
%% 主方法
[ KD,KM ] = GIPSim(interaction,1,1 );
Sk2 = beta * interaction' + beta^2 * (KM*interaction' + interaction'* KD);
F=Sk2';
F_ori(index_2,:) = F(index_2,:);
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
% [AUC1,AUPR1,Acc,Sen,Spe,Pre]=ROCcompute(F_ori,interaction,1);
 auc_ave = mean(auc);
 auc_std = std(auc);
 y(m,n)=auc_ave;
 n=n+1;
%  end 
 m=m+1;
 n=1;
%  end 
 o=toc;
 o
