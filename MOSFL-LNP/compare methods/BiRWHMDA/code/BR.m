function [Rt]=BR(A,Wrr,Wdd,L1,L2,alpha)

normWrr = normFun1(Wrr);
normWdd = normFun1(Wdd);

R0=A/sum(A(:));
Rt=R0;
nRtleft=0;
nRtright=0;
%bi-random walk on the heterogeneous network
for t=1:max(L1,L2)
    
    ftl = 0;
    ftr = 0;
    
    %random walk on the lncRNA similarity network
    if(t<=L1)
        nRtleft =(1- alpha) * normWrr * Rt + alpha*R0;
        ftl = 1;
    end
    %random walk on the disease similarity network
    if(t<=L2)
        nRtright = (1-alpha) *  Rt * normWdd + alpha*R0;
        ftr = 1;
    end
    
    %Rt: predictive association scores between each lncRNA-disease pair
    Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
    
end

end

