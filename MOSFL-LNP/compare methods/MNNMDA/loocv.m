clc;clear;
%% 原始处理
interaction = importdata('microbedata\interaction_disbiome.mat');
disSim = importdata('microbedata\disease_features1.txt');   %lncRNA 表达相似性
microSim = importdata('microbedata\microbe_features1.txt');  % disease 语义相似性
interaction_ori = interaction; 
n=1;
m=1;
%% 算法部分
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
maxiter = 300;
alpha = 1;
beta = 1;
tol1 = 2*1e-3;
tol2 = 1*1e-5;  
[dn,dr] = size(interaction);
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
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
    tic 
     u
    interaction(index(u))=0;
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
    F_ori(index(u))=F(index(u));
    interaction = interaction_ori;   
    toc
end
%% 结果
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');
figure;
[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );
[AUC,AUPR,Acc,Sen,Spe,Pre]=ROCcompute(F_ori,interaction,1); 
aupr=pr_cure(pre_label_score,label_y,'red');
% x(m,n)=aupr;
% y(m,n)=auc;
n=n+1;
% end
n=1;
m=m+1;
% end  
%end  




