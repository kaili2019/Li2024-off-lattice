clear
clc
% close all

exp_area =  573921; % from the average of colonies

% run simulations upto 5% of the experimental colony size
lwr_area = exp_area*0.95;
upr_area = exp_area*1.05;

% Total number of nutrients available 
N = 45000;

% parameters infered from mean posteior of ABC results

load("data/theta50_21100.mat")
mean50 = mean(theta50);

Telong = mean50(1);
p2sProb = mean50(2);
s2pProb = mean50(3);
pc = mean50(4);
pa = mean50(5);

%%
Ix = 1600;
Iy = 1200;

%%
for ii = 1:10
    [I,file_name] = run_off_lattice_v2(N,Telong,p2sProb,s2pProb,pc,pa,upr_area, ...
        lwr_area,exp_area,Ix,Iy);

    imwrite(I<0.5,"7_03Jun2024/AWRI796_50um_"+ii+".png")
end
