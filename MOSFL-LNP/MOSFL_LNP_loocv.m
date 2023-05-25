clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction.mat');
disSim = importdata('microbedata\disease_features.txt');   %disease 功能相似性
microSim = importdata('microbedata\microbe_features.txt');  % microbe 功能相似性
[nd,nm] = size(interaction);
interaction_ori = interaction; 
n=1;
m=1;
%% HOASF算法部分
% for alpha=0:0.5:4 
%   for lamda=[0,0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]     
alpha=2;
beta=1;
lamda =0.1;
Kd=32;
Km=5;
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
%% LNP主方法
F_ori = LNP(interaction_ori,KD,KM,nd,nm,lamda);
F_ori0 = F_ori;
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
    tic
     alpha;
     Km;
     lamda
     u
     
    interaction(index(u))=0;  
% computing liner neighborhood interaction profile kernel of lncRNAs and disease
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
F_ori(index(u))=F(index(u));
interaction = interaction_ori;
toc
end
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );
aupr=pr_cure(pre_label_score,label_y,'red');
% x(m,n)=aupr;
% y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  





