clc;clear;
tic
% format long
%% 原始数据
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %lncRNA 表达相似性
microSim = importdata('microbedata\microbe_features.txt');  % disease 语义相似性
[nd,nm] = size(interaction);
interaction_ori = interaction;
n = 1;
m=1;
c=1;
%% 数据处理方法
%  for alpha=0:0.2:4
%      for beta=0:0.2:4
% for lamda=0.1:0.1:1
%  for Kd=2:2:38
%      for Km=5:5:100     
% for lamda=0:0.1:1
alpha=2;
beta =1;      
lamda =0.1;
Kd = 2;
Km =10;
kd = Label_Propagation(interaction,0,Kd,'regulation2');
km = Label_Propagation(interaction',0,Km,'regulation2');
kd1  = cosSim( interaction); % disease cosine
km1  = cosSim( interaction');   % microbe cosine
w{1}=km;
w{2}=microSim;
w{3}=km1;
w1{1}=kd;
w1{2}=disSim;
w1{3}=kd1;
KM=MOSFL(w,Km,alpha,beta);
KD=MOSFL(w1,Kd,alpha,beta);

%% 主方法
   F_ori = LNP(interaction,KD,KM,nd,nm,lamda);
 %%
 F_ori_ori= F_ori;
%%  5-fold 交叉验证
auc = zeros(1,20);
% aupr = zeros(1,10);
% rng(0);
for k = 1:20
    tic
    alpha;
    beta;
    Kd
    Km
    lamda;
    k
    indices = crossvalind('Kfold', nd, 5 );
    interaction = interaction_ori;
    F_ori=F_ori_ori;
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(index_2,:) = 0;
       %%%计算得分矩阵
%%  数据处理
kd = Label_Propagation(interaction,0,Kd,'regulation2');
km = Label_Propagation(interaction',0,Km,'regulation2');
kd1  = cosSim( interaction); % disease cosine
km1  = cosSim( interaction');   % microbe cosine
w{1}=km;
w{3}=km1;
w1{1}=kd;
w1{3}=kd1;
KM=MOSFL(w,Km,alpha,beta);
KD=MOSFL(w1,Kd,alpha,beta);

%% 主方法
   F = LNP(interaction,KD,KM,nd,nm,lamda);
%%
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
 m=m+1;
 n=1;
%  end 
%  y1{c}=y;
%  c=c+1;
% end
toc
