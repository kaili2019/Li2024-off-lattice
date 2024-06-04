% this script make visualisation of the trace and histogram of results from  
% the ABC-MCMC algorithm based off the off-lattice model.
%
% Kai Li
% 26 May 2024

close all
clear
clc

% load data

load("data/theta50_21100.mat")
load("data/theta500_21400.mat")
load("data/thetaSW_22000.mat")

% plotting mean as marker 
figure('Position',[500,500,600,400])


A = [5, 2, 2, 2, 2]; 
B = [2, 2, 2, 5, 5];
param = ["$n^*$","$p_{ps}$","$p_{sp}$","$\gamma$","$p_a$"];

x = 0:0.01:1;

kde_50 = [0.06, 0.08, 0.08, 0.05, 0.03];
kde_500 = [0.03, 0.08, 0.05, 0.05, 0.03];
kde_SW = [0.05, 0.08, 0.05, 0.05, 0.03];

titles = ["(a)", "(b)", "(c)", "(d)", "(e)"];

meanWidth = 1;
distWidth = 1.5;

t = tiledlayout('flow','TileSpacing','compact');
for ii = 1:5
    % subplot(2,3,ii)
    nexttile
    hold on
    
    [f,xi] = ksdensity(theta50(:,ii),'Bandwidth',kde_50(ii)); 
    h1 = plot(xi,f,'color',[27,158,119]/255,'LineWidth',distWidth);
    t1 = plot(mean(theta50(:,ii),'omitnan'),0,'Marker','o','color',[27,158,119]/255);

    [f,xi] = ksdensity(theta500(:,ii),'Bandwidth',kde_500(ii)); 
    h2 = plot(xi,f,'color',[217,95,2]/255,'LineWidth',distWidth);
    t2 = plot(mean(theta500(:,ii),'omitnan'),0,'o', ...
        'Color',[217,95,2]/255);

    [f,xi] = ksdensity(thetaSW(:,ii),'Bandwidth',kde_SW(ii)); 
    h3 = plot(xi,f,'color',[117,112,179]/255,'LineWidth',distWidth);
    t3 = plot(mean(thetaSW(:,ii),'omitnan'),0, ...
        'o','Color',[117,112,179]/255);

    % priors
    p = plot(x,betapdf(x,A(ii),B(ii)),'Color',[.7 .7 .7],'LineStyle','--',...
             'LineWidth',distWidth);

    title(param(ii),'FontSize',24)
    xlabel(titles(ii),'FontSize',16)
    ylim([0,8])
    xlim([0,1])
    box on
end

hL = legend([h1,h2,h3,p,t3], ["AWRI 50$\mu M$", "AWRI 500$\mu M$",...
    "SW 50$\mu M$", "Prior","Mean"], ...
    'FontSize',16,'Box','off');
hL.Position(1:2) = [0.67,0.15];

% exportgraphics(gcf,"Fig9_histogram_v3.pdf")