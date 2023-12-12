clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %lncRNA 表达相似性
microSim = importdata('microbedata\microbe_features.txt');  % disease 语义相似性
interaction_ori = interaction; 
Z=zeros(1,450);
n=1;
m=1;
%% CSMF算法部分
lambdaM = 1;
lambdaD = 1;
lw = 0.7;
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
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
    tic 
     u
     
    interaction(index(u))=0;  
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
    F_ori(index(u))=F(index(u));
    Z(u)=F(index(u));
    interaction = interaction_ori;   
    toc
end
%% 结果
B=zeros(1,450);
for i=1:450
    if Z(i)>=0.3
        B(i)=1;
    end
end
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
% figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );
aupr=pr_cure(pre_label_score,label_y,'red');
% x(m,n)=aupr;
% y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  
%end  




