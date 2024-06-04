clear
clc
% close all

exp_area = 825329; % from the average of colonies

% run simulations upto 5% of the experimental colony size
lwr_area = exp_area*0.95;
upr_area = exp_area*1.05;

% Total number of nutrients available 
N = 65000;

% parameters infered from mean posteior of ABC results
load("data/thetaSW_22000.mat")
meanSW = mean(thetaSW);

Telong = meanSW(1);
p2sProb = meanSW(2);
s2pProb = meanSW(3);
pc = meanSW(4);
pa = meanSW(5);

%%
Ix = 1602;
Iy = 1418;

for ii = 1:10
    [I,file_name] = run_off_lattice_v2(N,Telong,p2sProb,s2pProb,pc,pa,upr_area, ...
        lwr_area,exp_area,Ix,Iy);
    imwrite(I<0.5,"7_03Jun2024/SW_"+ii+".png")
end
