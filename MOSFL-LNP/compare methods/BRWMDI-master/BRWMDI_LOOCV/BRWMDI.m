function [F] = BRWMDI(Y, dsymsim, K, Nt, alpha, Il, Ir)
    %GIP parameters
    gamadm = 1;
    gamadd = 1 ;
    g1=-9;
    g2=log(9999);
    [nd,nm] = size(Y);
    %calculate gamad for Gaussian kernel calculation
    for i=1:nd
        sd(i)=norm(Y(i,:))^2;
    end
    gamad=nd/sum(sd')*gamadd; 
    %calculate gamam for Gaussian kernel calculation
    for i=1:nm
        sm(i)=norm(Y(:,i))^2;
    end
    gamam=nm/sum(sm')*gamadm; 
    %calculate Gaussian kernel for the similarity between disease: kd
    for i=1:nd
        for j=1:nd
            kd(i,j)=exp(-gamad*(norm(Y(i,:)-Y(j,:)))^2);
        end
    end   
    %calculate Gaussian kernel for the similarity between microbe: km
    for i=1:nm
        for j=1:nm
            km(i,j)=exp(-gamam*(norm(Y(:,i)-Y(:,j)))^2);
        end
    end 
        
    
    %Similarity network fusion
%     microSim = importdata('microbedata\microbe_features.txt');  % microbe 功能相似性
    kd = SNF({kd,dsymsim}, K, Nt);
%     kd = Label_Propagation(Y,0,32,'regulation2');
%     km = Label_Propagation(Y',0,5,'regulation2');
%     kd1  = cosSim( Y); % disease cosine
%     km1  = cosSim( Y');   % microbe cosine
%     Md = HSCPES({kd,dsymsim,kd1},32,2,1);
%     Mm = HSCPES({km,microSim,km1},5,2,1);
    [km] = logisticfunction(km,g1,g2);
    Md = kd;
    Mm = km;
    
    Y = Y / sum(Y(:));
    R = Y;
    Max_Iter = max([Ir, Il]);
    mnumber = size(Mm,1);
    dnumber = size(Md,1);
    F = zeros(dnumber,mnumber);

    for i = 1 : size(Mm,2)
        sumcoltmp = sum(Mm(:, i));
        if sumcoltmp > 0
            Mm(:, i) = Mm(:, i) / sumcoltmp;
        end
    end
    for i = 1 : size(Md, 2)
        sumcoltmp = sum(Md(:, i));
        if sumcoltmp > 0
            Md(:, i) = Md(:, i) / sumcoltmp;
        end
    end


    for i = 1 : Max_Iter
        rflag = 0;
        lflag = 0;

        if i <= Il
            lflag = 1;
            LR = (1 - alpha) * Md * R + alpha * Y; 
        end

        if i <= Ir
            rflag = 1;
            RR = (1 - alpha) * R * Mm + alpha * Y;
        end

        R = (rflag * RR + lflag * LR) / (rflag + lflag); 
    end


    for i = 1 : dnumber
        for j = 1 : mnumber
            F(i,j)=R(i,j);
        end
    end