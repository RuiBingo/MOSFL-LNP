%%  ����˵��    20181215  version1@lotus
%   aupr =pr_cure��Ԥ�����,ԭʼ��׼��,colour)
%   ����PR���ߺ�PR���������aupr
%   PR���ߣ�precison-recall curve��

function aupr =pr_cure(output,original,colour)

    %%  ��������1
    %{
        colour='red';
        %label_y=[1 1 1 0 1 0 1 1 0 1];
       % deci=[0.89 0.7 0.87 0.32 0.50 0.14 0.44 0.59 0.74 0.99];
        label_y=[0 0 0 0 1 1 1 1 1 1 ];
        deci=[0.99 0.2 0.3 0.4  0.5 0.6 0.7 0.8 0.9 0.95];
    %}

    %%  ��������2
    %{
    colour='red';
    original=textread('pro_lnc_matrix.txt');
    original=original(:);
    load('lncpro_LPBNI_LOOCV.mat'); 
    output=auc_matrix;
    output=output(:);
    %}

    %% ��Ԥ��������deci�������򣬱�׼��roc
    [threshold,ind] = sort(output,'descend');  %[��ֵ���±�]����Ԥ�������������
    roc_y = original(ind);    %����ֵԤ������Ӧ�ı�׼��

    %% ��x��recall�ĸ����㣬��y��precison�ĸ����㡣��PR���������aupr
    P=[1:length(roc_y)]';   %ʵ��������(TP+FP)��������Ԥ��Ϊ���ĸ�������Ϊ��ֵ���ǽ���������ֵ��Ӧ���±꼴��TP+FP��
    stack_x = cumsum(roc_y == 1)/sum(roc_y == 1); %x�᣺TPR=recall=TP/(TP+FN)=Ԥ��Ϊ��������/��������
    stack_y = cumsum(roc_y == 1)./P; %y�᣺precision=TP/(TP+FP)=Ԥ��Ϊ��������/����Ԥ��Ϊ��
    %stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
    %stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
    aupr=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));  %PR���������

    %% ��PR����
    % subplot(2,2,1);   %�ѻ�ͼ���ڷֳ����������Ŀ�����Ȼ����ÿ������ֱ���ͼ���ڵ�һ���ͼ
    %figure;
    plot(stack_x,stack_y,colour);
    xlabel('recall');
    ylabel('precision');
    title(['ROC curve of (AUC = ' num2str(aupr) ' )']);
end



%% 
%����myACC_1��ACC��F1_score�ȵ����ַ���ʾ��
    %myACC_1( deci,label_y,'sp0.99' ); %  method��'sp0.99'��'sp0.95'�� 'youden'��'accMAX'����
    %myACC_1( deci,label_y,'sp0.95' );
    %myACC_1( deci,label_y,'youden' );
    %myACC_1( deci,label_y,'accMAX' );
