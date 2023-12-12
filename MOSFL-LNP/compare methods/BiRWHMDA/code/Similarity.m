function [kl,kd]=Similarity(interaction,gamall,gamadd)
% Similarity(1,1)
%construct the adjacency matrix for the microbe-disease association network and calculate microbe similarity and disease similarity
  %A: Binary relations between diseases and microbes, 1st column:disease, 2nd column:microbe
%   clc,clear;
%   A=textread('knowndiseasemicrobeinteraction.txt');
% gamall=1;
% gamadd=1;
  % nd:the number of diseases
  % nl:the number of microbes
  % pp:the number of known microbe-disease associations
%   nd=max(A(:,1)); 
%   nl=max(A(:,2));
%   [pp,qq]=size(A);

  %interaction: adjacency matrix for the microbe-disease association network
  %interaction(i,j)=1 means there is a known association between microbe j and disease i
%   for i=1:pp
%         interaction(A(i,1),A(i,2))=1;
%   end
  
%   save ('interaction.mat','interaction');
  [nd,nl] = size(interaction);
  %calculate gamal for Gaussian kernel calculation
  for i=1:nl
      sl(i)=norm(interaction(:,i))^2;
  end
  gamal=nl/sum(sl')*gamall;

  %calculate gamad for Gaussian kernel calculation
  for i=1:nd
        sd(i)=norm(interaction(i,:))^2;
  end
  gamad=nd/sum(sd')*gamadd;
    
  %calculate Gaussian interaction profile kernel similarity for microbes: kl
  for i=1:nl
       for j=1:nl
           kl(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
       end
  end 
  
%   save ('microbesimilarity.mat','kl');
  
  %calculate Gaussian interaction profile kernel similarity for disease: pkd
  for i=1:nd
      for j=1:nd
          pkd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
      end
  end
    
  %logistic function transformation for disease similarity
  for i=1:nd
      for j=1:nd
          kd(i,j)=1/(1+exp(-15*pkd(i,j)+log(9999)));
      end
  end
  
%   save ('diseasesimilarity.mat','kd');

end



        
        
        
    
   



