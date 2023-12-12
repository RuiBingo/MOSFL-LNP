clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %lncRNA 表达相似性
microSim = importdata('microbedata\microbe_features.txt');  % disease 语义相似性
interaction_ori = interaction; 
Z=zeros(1,450);
n=1;
m=1;
%% 算法部分
    g1=-15;
    g2=log(9999);
[kl,kd]=Similarity(interaction,1,1);
% [kd] = logisticfunction(kd,g1,g2);
% Rt_ori=BiRWHMDA(interaction,KD,KM,0.4,2,2);
Rt_ori=BiRWHMDA(interaction,kd,kl,0.4,2,2);
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
    tic 
     u
     
    interaction(index(u))=0;  
    [kl,kd]=Similarity(interaction,1,1);
    Rt=BiRWHMDA(interaction,kd,kl,0.4,2,2);
    Rt_ori(index(u))=Rt(index(u));
    Z(u)=Rt(index(u));
    interaction = interaction_ori;   
    toc
end
%% 结果
B=zeros(1,450);
for i=1:450
    if Z(i)>=10^(-5)
        B(i)=1;
    end
end
pre_label_score = Rt_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
% figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( Rt_ori(:),interaction_ori(:),'sp0.99' );
aupr=pr_cure(pre_label_score,label_y,'red');
% x(m,n)=aupr;
y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  
%end  




