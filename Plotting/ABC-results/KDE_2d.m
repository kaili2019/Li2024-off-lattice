% this script make visualisation of histogram of the ABC-MCMC of the 
% off-lattice model. KDE applies the kenel density estimator over the 
% accpeted samples producing n* vs gamma. 
% 
% Kai Li
% 4 June 2024

close all
clear
clc

% load data

load("data/theta50_21100.mat")
load("data/theta500_21400.mat")
load("data/thetaSW_22000.mat")


%% plotting
 
% Define grid
xgrid=linspace(0,1,100);
ygrid=linspace(0,1,100);
[x1,y1] = meshgrid(xgrid, ygrid);

% Perform kernel density estimate
xi = [x1(:) y1(:)];
close all
df{1} = theta50;
df{2} = theta500;
df{3} = thetaSW;



for ii = 1:3
    figure('Position',[400,500,400,300])
    x = df{ii}(:,1);
    y = df{ii}(:,4);
    
    [f,ep]=ksdensity([x y],xi,'Kernel','box','Bandwidth',0.05); % remove the outputs to see a 3D plot of the distribution
    
    % format data in matrix for contourf and plot
    X = reshape(ep(:,1),length(xgrid),length(ygrid));
    Y = reshape(ep(:,2),length(xgrid),length(ygrid));
    Z = reshape(f,length(xgrid),length(ygrid));
    contourf(X,Y,Z)

    xlabel("$n^*$",'FontSize',36)
    ylabel("$\gamma$",'FontSize',36)
    hold on
    plot(mean(x,'omitnan'),mean(y,'omitnan'),'Marker','square','MarkerFaceColor','white','MarkerEdgeColor','k',...
        'MarkerSize',12)
    axis equal
    colorbar
    % axis square
    % scatter(x,y,'r*')
    axis equal

    % save figure
%     if ii == 1
%         % title("AWR 796 50 $\mu M$",'FontSize',24)
%         exportgraphics(gcf,'ABC_kde50um3.pdf','ContentType','vector')
%     elseif ii == 2
%         % title("AWR 796 500 $\mu M$",'FontSize',24)
%         exportgraphics(gcf,'ABC_kde500um3.pdf','ContentType','vector')
%     else
%         % title("Simi White 50 $\mu M$",'FontSize',24)
%         exportgraphics(gcf,'ABC_kdeSW3.pdf','ContentType','vector')
%     end
end




