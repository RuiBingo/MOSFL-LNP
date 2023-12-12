clc;clear;
tic

%% 原始数据
interaction = importdata('microbedata\interaction.mat');
[nl,nd] = size(interaction);
interaction_ori = interaction;
n = 1;
m=1;
%% 主方法
[kl,kd]=Similarity(interaction,1,1);
F_ori=BiRWHMDA(interaction,kd,kl,0.4,2,2);
index=find(interaction_ori==1);
%%  5-fold 交叉验证
auc = zeros(1,20);
% rng(0);
for k = 1:20
    k
    indices = crossvalind('Kfold', nd, 5 );
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(:,index_2) = 0;
       %%%计算得分矩阵
%% 主方法
[kl,kd]=Similarity(interaction,1,1);
F=BiRWHMDA(interaction,kd,kl,0.4,2,2);
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
%  end 
 m=m+1;
 n=1;
%  end 
 o=toc;
 o
