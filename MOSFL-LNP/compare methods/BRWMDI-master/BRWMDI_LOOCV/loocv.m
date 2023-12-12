clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction.mat');
% disSim = importdata('dsymsim.mat'); % disease 症状相似性
disSim = importdata('microbedata\disease_features.txt');
interaction_ori = interaction; 
n=1;
m=1;
%% 算法部分
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
K = 4;%number of neighbors
Nt = 5; %Number of Iterations       
%BRWH parameters, alpha:decay facor; 
%Il:maximum iteration numbers of disease network
%Ir:maximum iteration numbers of microbe network
alpha = 0.5;
Il = 1;
Ir = 2;
F_ori = BRWMDI(interaction, disSim, K, Nt, alpha, Il, Ir);
% allresult1(microbe,disease,interaction_ori',F_ori')
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
    tic 
     u
    interaction(index(u))=0;
    F = BRWMDI(interaction, disSim, K, Nt, alpha, Il, Ir);
    F_ori(index(u))=F(index(u));
    interaction = interaction_ori;   
    toc
end
%% 结果
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
% figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );
% [auc,AUPR,Acc,Sen,Spe,Pre]=ROCcompute(F_ori,interaction_ori,1); 
aupr=pr_cure(pre_label_score,label_y,'red');
% x(m,n)=aupr;
y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  
%end  




