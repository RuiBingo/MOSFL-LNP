clc;clear;
tic
% format long
%% 原始数据
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %lncRNA 表达相似性
microSim = importdata('microbedata\microbe_features.txt');  % disease 语义相似性
[dn,dr] = size(interaction);
interaction_ori = interaction;
n = 1;
m=1;
%% 主方法
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
maxiter = 300;
% for alpha = [0.1,1,10,100]
%     for beta = [0.1,1,10,100,1000]
alpha = 1;
beta = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;  
[ GD,GM ] = GIPSim(interaction,1,1 );
    HD=zeros(nd,nd);
    HM=zeros(nc,nc);
    for i=1:nd
        for j=1:nd
            if disSim(i,j)>0
                HD(i,j)=(disSim(i,j)+GD(i,j))/2;
            else
                HD(i,j)=GD(i,j);
            end
        end
    end
    for i=1:nc
        for j=1:nc
            if microSim(i,j)>0
                HM(i,j)=(microSim(i,j)+GM(i,j))/2;
            else
                HM(i,j)=GM(i,j);
            end
        end
    end
    T = [HM, interaction'; interaction, HD];
    [t1, t2] = size(T);
    trIndex = double(T ~= 0);
    [WW,iter] = DRMNN(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
    F_ori = WW((t1-dn+1) : t1, 1 : dr);
% index=find(interaction_ori==1);
%%  5-fold 交叉验证
auc = zeros(1,20);
% aupr = zeros(1,100);
% rng(0);
for k = 1:20
    alpha;
    beta;
    k
    indices = crossvalind('Kfold', nc, 5 );
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(:,index_2) = 0;
       %%%计算得分矩阵
%% 主方法
[ GD,GM ] = GIPSim(interaction,1,1 );
    HD=zeros(nd,nd);
    HM=zeros(nc,nc);
    for i=1:nd
        for j=1:nd
            if disSim(i,j)>0
                HD(i,j)=(disSim(i,j)+GD(i,j))/2;
            else
                HD(i,j)=GD(i,j);
            end
        end
    end
    for i=1:nc
        for j=1:nc
            if microSim(i,j)>0
                HM(i,j)=(microSim(i,j)+GM(i,j))/2;
            else
                HM(i,j)=GM(i,j);
            end
        end
    end
    T = [HM, interaction'; interaction, HD];
    [t1, t2] = size(T);
    trIndex = double(T ~= 0);
    [WW,iter] = DRMNN(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
    F = WW((t1-dn+1) : t1, 1 : dr);
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
% [AUC1,AUPR1,Acc,Sen,Spe,Pre]=ROCcompute(F_ori,interaction,1);
 auc_ave = mean(auc);
 auc_std = std(auc);
 y(m,n)=auc_ave;
 n=n+1;
% end 
 m=m+1;
 n=1;
% end 
toc

