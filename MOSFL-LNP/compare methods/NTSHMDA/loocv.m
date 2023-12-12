clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction.mat');
% disSim = importdata('microbedata\disease_features1.txt');   %lncRNA 表达相似性
% microSim = importdata('microbedata\microbe_features1.txt');  % disease 语义相似性
interaction_ori = interaction; 
Z=zeros(1,450);
n=1;
m=1;
%% CSMF算法部分
gamma=0.7;
% for beta2=0.1:0.1:0.9
phi=0.9;
delte=0.3;
beta1=0.4;
beta2=0.2;
% beta2=1-beta1;
% F_ori=NTSHMDA(interaction,disSim,microSim,gamma,phi,delte,beta1,beta2);
F_ori=NTSHMDA(interaction,gamma,phi,delte,beta1,beta2);
% allresult1(microbe,disease,interaction_ori',F_ori')
%%  留一交叉验证（LOOCV
index=find(interaction_ori==1); 

for u=1:length(index)
    tic 
     u
%      beta2
    interaction(index(u))=0;  
    F=NTSHMDA(interaction,gamma,phi,delte,beta1,beta2);
%     F=NTSHMDA(interaction,disSim,microSim,gamma,phi,delte,beta1,beta2);
    F_ori(index(u))=F(index(u));
    Z(u)=F(index(u));
    interaction = interaction_ori;   
    toc
end
B=zeros(1,450);
for i=1:450
    if Z(i)>=10^(-6)
        B(i)=1;
    end
end
%% 结果
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
% figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );
aupr=pr_cure(pre_label_score,label_y,'red');
x(m,n)=aupr;
y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  
%end  




