clc;
clear;
BRWMDI_LOOCV();
overallauc=positiontooverallauc();
save overallauc overallauc
c=mean(overallauc)
d=std(overallauc)