% this script make visualisation of the trace plots from the ABC-MCMC 
% algorithm based off the off-lattice model. T
% 
% Kai Li
% 04 June 2024

close all
clear
clc

% load data

load("data/theta50_21100.mat")
load("data/theta500_21400.mat")
load("data/thetaSW_22000.mat")

%% traceplots

figure('Position',[0,200,1300,300])
plot(thetaSW,'LineWidth',1)
legend({'$n^*$','$P_{ps}$','$P_{sp}$','$\gamma$','$P_a$'},'FontSize',16)
xlabel("Interation $i$",'FontSize',24)
ylabel("$\theta$",'FontSize',24)
xlim([0,20000])
% title("Simi White $50 \mu M$",'FontSize',24)
% exportgraphics(gcf,'SW_traceplot.pdf','ContentType','vector')

%%
figure('Position',[0,200,1300,300])
plot(theta50,'LineWidth',1)
legend({'$n^*$','$P_{ps}$','$P_{sp}$','$\gamma$','$P_a$'},'FontSize',16)
xlabel("Interation $i$",'FontSize',24)
ylabel("$\theta$",'FontSize',24)
xlim([0,20000])
% title("Simi White $50 \mu M$",'FontSize',24)
% exportgraphics(gcf,'50um_traceplot.pdf','ContentType','vector')

%%
figure('Position',[0,200,1300,300])
plot(theta500,'LineWidth',1)
legend({'$n^*$','$P_{ps}$','$P_{sp}$','$\gamma$','$P_a$'},'FontSize',16)
xlabel("Interation $i$",'FontSize',24)
ylabel("$\theta$",'FontSize',24)
xlim([0,20000])
% title("Simi White $50 \mu M$",'FontSize',24)
% exportgraphics(gcf,'500um_traceplot.pdf','ContentType','vector')

%%
clc
disp("50um ESS")
multiESS(theta50)

disp("500um ESS")
multiESS(theta500)

disp("SW ESS")
multiESS(thetaSW)