Using computer technology to study bioinformatics. Code and Datasets of Paper--Predicting microbe-disease associations based on a linear neighborhood label propagation method with multi-order similarity fusion learning

compare methods：The compare methods-file contains all the code for the comparison algorithms

Input:

microbe_features.txt (microbe similarity metrix)

interaction.mat (microbe-disease associations)

disease_features.txt (disease similarity metrix)

Output:

Predictive scoring matrix

Main code description：

MOSFL.m -- multi-order similarity fusion learning

Label_Propagation.m -- Calculating linear neighborhood similarity

LNP.m -- linear neighborhood label propagation

MOSFL_LNP_loocv.m -- leave-one-out cross-validation frameworks

MOSFL_LNP_5_fold.m -- 5-fold validation frameworks

cv_2.m -- 5-fold CVS2

cv_3.m -- 5-fold CVS3
