function F = LNP(interaction,SD,SM,nd,nm,alpha)

LPD = (1-alpha)*pinv(eye(size(SD,1))-alpha*SD)*interaction;
LPM = (1-alpha)*pinv(eye(size(SM,1))-alpha*SM)*interaction';
LPM = LPM';

% 最值归一化
for i=1:nd

     MAX_D = max(LPD(i,:));
    MIN_D = min(LPD(i,:));

        LPD(i,:) = (LPD(i,:) - MIN_D) / (MAX_D - MIN_D); 
    

end

 for i =1:nm
     MAX_M = max(LPM(:,i));
     MIN_M = min(LPM(:,i));
 
     LPM(:,i) = (LPM(:,i) - MIN_M) / (MAX_M - MIN_M);
  
 end

beta = 1/2;
F = beta * LPM + (1-beta) * LPD;

