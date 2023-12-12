clc;clear;
tic
% format long
%% 原始数据
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %疾病 相似性

[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
interaction_ori = interaction;
 n = 1;
 m=1;
%% 数据处理方法
%  for K=10:10:100
%      for K3=10:10:100 
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
K = 4;%number of neighbors
Nt = 5; %Number of Iterations       
%BRWH parameters, alpha:decay facor; 
%Il:maximum iteration numbers of disease network
%Ir:maximum iteration numbers of microbe network
alpha = 0.9;
% for alpha = 0.1:0.1:0.9
%     for Ir = 1:6
Il = 1;
Ir = 2;
F_ori = BRWMDI(interaction, disSim, K, Nt, alpha, Il, Ir);
 %%
 F_ori_ori= F_ori;
 index=find(interaction_ori==1);
%%  5-fold 交叉验证
auc = zeros(1,100);
% rng(0);
for k = 1:100
    tic
    k
    indices = crossvalind('Kfold', length(index), 5 );
    interaction = interaction_ori;
    F_ori=F_ori_ori;
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(index(index_2)) = 0;
       %%%计算得分矩阵
%% 主方法
   F = BRWMDI(interaction, disSim, K, Nt, alpha, Il, Ir);
%%
      F_ori(index(index_2)) = F(index(index_2));
      interaction = interaction_ori;
end
%% 画auc曲线 
    pre_label_score = F_ori(:);
    label_y = interaction_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');
%     [ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.95' );
%     figure;
%     aupr(k)=pr_cure(pre_label_score,label_y,'red');
    toc
end
%% 
 auc_ave = mean(auc);
 auc_std = std(auc);
% x(m,n)=aupr;
%   x(m,n)=gama;
   y(m,n)=auc_ave;
   n=n+1;
%  end 
%  m=m+1;
%  n=1;
%  end 
toc
 
