%%  ����˵��    20181215  version1@lotus
%  [ACC,PRE,SEN,F1_score,MCC] = myACC_1(Ԥ�����,ԭʼ��׼��,method)
%  method='sp0.99'or'sp0.95'or 'youden'or 'accMAX'
        %'sp0.99'or'sp0.95'����sp=0.99��0.95��Ϊ��ֵ��� ��%�ҵ���ӽ�FPR_value=0.01��0.05��λ��
        %'youden'����YoudenԼ��ָ������ֵ,����TPR-FPRȡ���ֵʱ��Լ������������Ͻǵĵ�ʱ����ֵ
        %'accMAX'����ACC���ֵ����ֵ

function [ACC,PRE,SEN,F1_score,MCC] = myACC_1( output,original,method )
      
    %{
    %��������  
      original=textread('pro_lnc_matrix.txt');
      load('lncpro_LPBNI_LOOCV.mat'); 
      output=auc_matrix;
     %}
   
    %% ����Ƿ��ѽ�������һ��
    output=output(:);
    original=original(:);

    %% ͳ������͸������
    count1=length(find(original==1));  %�������=TP+FN=����+����
    count0=length(find(original==0));  %�������=FP+TN=����+����

    [tpr,fpr,thresholds] = roc(original(:)',output(:)');  %[�����ʣ�������,��ֵ]
    TP=count1*tpr;  %TP�������������
    FP=count0*fpr;  %FP�����������
    ACC= (TP+count0-FP)/(count1+count0);  %׼ȷ��ACC=(TP+TN)/all=(����+����)/all  ,����TN=count0-����FP

    %% ѡ��һ�ַ���������ֵ  
        % method='sp0.99'��'sp0.95'�� 'youden'��'accMAX'����
     
        %�趨sp��ֵ��ACC
        if method=='sp0.99' 
            sp_value=0.99 ;  %  sp=0.99��1-sp=FPR=0.01��
            [~,position]=min(abs( fpr-(1-sp_value)));  %�ҵ���ӽ�FPR_value��λ��

        %�趨sp��ֵ��ACC�� 
        elseif method=='sp0.95'
            sp_value=0.95 ;   % sp=0.95��1-sp=FPR=0.05��
            [~,position]=min(abs( fpr-(1-sp_value)));  %�ҵ���ӽ�FPR_value��λ��

        %��YoudenԼ��ָ������ֵ
        elseif method=='youden'
            [~,position]=max(abs(tpr-fpr));  %��Լ��ָ������ֵ,����TPR-FPRȡ���ֵʱ��

        %��ACC���ֵ����ֵ
        elseif method=='accMAX'
            [~,position]=max(ACC); 
        end

        %��һ�㲻�ã����ض���ֵ��ACC
        %thre_value=0.05;  %�趨��ֵ��ֵ��ACC
        %[~,position]=min(abs(thresholds-thre_value));

    %% ��ǰ��ֵ��ACC�Ƚ��
    thre=thresholds(position); %��ǰ��ֵ
    TP=count1*tpr(position);  %TP�������������
    FP=count0*fpr(position);  %FP�������������
    FN=count1-TP;
    TN=count0-FP;

    %% ������PRE/SEN/ACC/F1_score/MCC
    PRE=TP/(TP+FP)  ;  %��ȷ��PRE=TP/TP+FP
    SEN=tpr(position); %���ж�SEM=������TPR        %����SEN=TP/(count1)
    ACC= (TP+count0-FP)/(count1+count0);  %׼ȷ��ACC=(TP+TN)/all=(����+����)/all  ,����TN=count0-����FP
    F1_score=2*TP/(count1+TP);  %F1����=2P*R/(P+R)����ȷ��P �����ж�R
    MCC=(TP*TN -FP*FN)/sqrt(count1*count0*(TP+FP)*(FN+TN));
    fprintf([ method,'��ACC=',num2str(ACC), '   PRE=',num2str(PRE), '   SEN=',num2str(SEN), '  F1 score=',num2str(F1_score), '   MCC=',num2str(MCC),  '\n']);
    
    %auc=roc_1(output(:),original(:),'red') %��һ��ROC
end

