function Rt=BiRWHMDA(A,Wdd,Wrr,alpha,l,r)
%BiRWHMDA(0.4,2,2)
%bi-random walk on the heterogeneous network to calculate association scores for each microbe-disease pair. 

  %Wrr: adjacency matrix of the microbe similarity network
  %Wdd: adjacency matrix of the disease similarity network
  %A: adjacency matrix of the microbe-disease association network
%   Wrr = importdata('microbesimilarity.mat');
%   Wdd = importdata('diseasesimilarity.mat');
%   A = importdata('interaction.mat');

  %normFun: Laplacian normalization for microbe similarity and disease similarity
  normWrr = normFun(Wrr);
  normWdd = normFun(Wdd);

  %R0: initial probability
  R0 = A/sum(A(:));
  Rt = R0;
  
  %bi-random walk on the heterogeneous network
  for t=1:max(l,r)
      
      ftl = 0;
      ftr = 0;
    
      %random walk on the microbe similarity network
      if(t<=l)
         nRtleft = alpha *  Rt * normWrr  + (1-alpha)*R0;
         ftl = 1;
      end   
      
      %random walk on the disease similarity network
      if(t<=r)
         nRtright = alpha * normWdd * Rt + (1-alpha)*R0;
         ftr = 1;
      end
    
      %Rt: predictive association scores between each microbe-disease pair
      Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
      
  end

%   save ('Result.mat','Rt');
  
end

