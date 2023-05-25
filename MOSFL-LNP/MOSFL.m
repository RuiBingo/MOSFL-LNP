function [W]=MOSFL(Wall,K,alpha,beta)

Wall1=Wall;
C = length(Wall);
[m,n]=size(Wall{1});
rho=0.5;

%% normalize
for i = 1 : C
    newW1{i} = Wall{i}./repmat(sum(Wall{i},1),n,1); 
    TF = isnan(newW1{i});
    newW1{i}(find(TF(:,:)==1)) = 0;       
end

for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end
%% initialize Omega and W 
for i = 1 : C
    Omega{i}=zeros(m,n);
end
Wsum=zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
P=Wsum/C;

ITER=0;
dual_norm=0;
%last=0;
    for i = 1 : C
        dual_norm = dual_norm+(norm(Wall{i}-P))/C;
    end
ez = 1;
 while(max(dual_norm/rho^2,ez)>0.000001)
    % update W1,...WC
    for i = 1 : C      
        %Wi-HiSj
        tmp2=zeros(n);
        for j=1:C
            tmp2=tmp2+(2*alpha*newW1{i}*newW{j}');
        end
        %Wi-Hj
        tmp1=zeros(n);
        for j=1:C
            tmp1=tmp1+2*newW1{j};
        end
        %Wi-SiKSj^T
        tmp=zeros(n);
        for j=1:C
           wtmp = FindDominateSet(Wsum,round(K));
           tmp=tmp+(2*beta*newW{i}*wtmp*newW{j}');
        end
        Wall{i}=(tmp1+tmp2+tmp+rho*P-rho*Omega{i})./((2*C+2*alpha*C+2*beta*C+rho)*ones(n));
    end

    % update W
    Wsum=zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
     end
    prev_P=P;
    P=Wsum/C;    
    %update Omega
    for i = 1 : C
        Omega{i}=Omega{i}+Wall{i}-P;
    end
    temp=0;
    for i = 1 : C
        temp = temp+(norm(Wall{i}-P))/C;
    end
    dual_norm(ITER+1) = temp;
    ez(ITER+1) = norm(prev_P-P);
    if dual_norm(ITER+1)>100*ez(ITER+1)
        rho = rho*0.8;
    elseif (ez(ITER+1)>.1*dual_norm(ITER+1) )
        rho = 10*rho;
    end
    ITER=ITER+1;
    if ITER>10000
        ITER
        break;
    end
end

W=P;
W=W./repmat(sum(W,2),1,n);
W = (W +W')/2;
end

function newW = FindDominateSet(W,K)
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
%I1 KNN Index
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
TF = isnan(newW);
newW(find(TF(:,:)==1)) = 0;
%newW=newW.*2;
clear IW1;
clear IW2;
clear temp;
end