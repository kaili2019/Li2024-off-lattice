clear
clc
% close all

exp_area =  1052348; % from the average of colonies

% run simulations upto 5% of the experimental colony size
lwr_area = exp_area*0.95;
upr_area = exp_area*1.05;

% Total number of nutrients available 
N = 90000;

% parameters infered from mean posteior of ABC results

load("data/theta500_18900.mat")
mean500 = mean(theta500);

Telong = mean500(1);
p2sProb = mean500(2);
s2pProb = mean500(3);
pc = mean500(4);
pa = mean500(5);

% simulate colonies
Ix = 1600;
Iy = 1500;

for ii = 7:10
    [I,file_name] = run_off_lattice_v2(N,Telong,p2sProb,s2pProb,pc,pa,upr_area, ...
        lwr_area,exp_area,Ix,Iy);
    imwrite(I<0.5,"7_03Jun2024/AWRI796_500um_"+ii+".png")
end
