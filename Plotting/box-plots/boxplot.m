% this script calculates the summary statistic of the final time of
% colonies under different conditions and produces boxplot using v2 of the
% summary statistic
% 
% Kai Li
% 21 April 2024

clear
clc
close all

% AWRI 796
AWRI_r = nan(3,14);
AWRI_f = nan(3,14);
AWRI_s = nan(3,14);

tic
kk = 1;
for cc = [1 5 7]

    if cc == 5 % 500um 
        ss = 2:14;
    else
        ss = 1:14;
    end

    for ii = ss   
        
        if cc == 1
            I = imread("Data/AWRI 796 "+50+"um/AWRI 796 "+50 ...
                +"um s"+ii+" 233h 4X binary.tif")<0.5;
        elseif cc == 5
            I = imread("Data/AWRI 796 "+500+"um/AWRI 796 "+500 ...
                +"um s"+ii+" 233h 4X binary.tif")<0.5;
        elseif cc == 7
            I = imread("Data/Simi White 50um/Simi White s"+ii+" 237h 4X binary.tif")<0.5;
        else
            disp("no colony detected")
        end
        
        AWRI_r(kk,ii) = get_norm_radius_ratio_v2(I);
        AWRI_f(kk,ii) = get_norm_filament_ratio_v2(I);
        AWRI_s(kk,ii) = get_norm_subbranch_v2(I);
    end
    kk=kk+1;
end
toc

%% option 1 visualisation

figure('Position',[300,300,400,300])
boxchart(AWRI_s')
hold on
plot(1:3,mean(AWRI_s,2,'omitnan'),'Marker','+','LineStyle','none')
ylim([0,1])
set(gca,'XTickLabel',{'AWRI $50\mu M$','AWRI $500\mu M$','SW $50\mu M$'},'YTick',0:0.2:1);
ylabel("$I_B$",'FontSize',36)
axis square
box on
% exportgraphics(gcf,'Fig8_ss_Ib_boxplots.pdf')

%%
pause(0.1)
figure('Position',[800,300,400,300])
boxchart(AWRI_f');
hold on
plot(1:3,mean(AWRI_f,2,'omitnan'),'Marker','+','LineStyle','none')
ylim([0,1])
set(gca,'XTickLabel',{'AWRI $50\mu M$','AWRI $500\mu M$','SW $50\mu M$'},'YTick',0:0.2:1);
ylabel("$I_F$",'FontSize',36)
axis square
box on
% exportgraphics(gcf,'Fig8_ss_If_boxplots.pdf')

%%
pause(0.1)
figure('Position',[1200,300,400,300])
boxchart(AWRI_r')
hold on
plot(1:3,mean(AWRI_r,2,'omitnan'),'Marker','+','LineStyle','none')
ylim([0,1])
set(gca,'XTickLabel',{'AWRI $50\mu M$','AWRI $500\mu M$','SW $50\mu M$'},'YTick',0:0.2:1);
ylabel("$I_R$",'FontSize',36)
axis square
box on

% exportgraphics(gcf,'Fig8_ss_Ir_boxplots.pdf')

%% calculate mean and sd of individual environment groups
clc
disp("branch")
disp([mean(AWRI_s,2,'omitnan'), std(AWRI_s,0,2,'omitnan')])

%%
clc
disp("filament")
disp([mean(AWRI_f,2,'omitnan'),std(AWRI_f,0,2,'omitnan')])

%%
clc
disp("radii")
disp([mean(AWRI_r,2,'omitnan'),std(AWRI_r,0,2,'omitnan')])