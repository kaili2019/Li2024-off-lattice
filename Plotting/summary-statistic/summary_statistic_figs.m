% this script makes the skeletonized image and cropped colony figures using
% CSR radius 
% 
% 18 April 2024
% Kai Li

clear 
clc
close all

% read in data
AWRI796_50um = imread("AWRI 796 50um s5 233h 4X binary.tif") < 0.5;
SW = imread("Simi White s14 237h 4X binary.tif") < 0.5;
AWRI796_500um = imread("AWRI 796 500um s14 233h 4X binary.tif") < 0.5;

colony = AWRI796_50um;
colony(end:end+100,:) = 0;

cut = 100; % padding binary figure

%% get radii
[x0,y0,rmax,rcsr] = get_csr_radius_v2(colony,4);

% get a matrix to crop out the centre of the colony based on csr radius
M = zeros(size(colony,1),size(colony,2));
M(round(y0),round(x0)) = 1;
R = bwdist(M);
C = R >= rcsr;

% crop a quarter of the colony
T = zeros(size(colony,1),size(colony,2));
T(1:round(y0),1:round(x0)) = 1;  

% skeletonize
skeleton = bwmorph(colony.*C, 'skel', Inf);

%% plotting skeletonized colony fig (b)

colony_crop = colony.*C;

if size(skeleton,2) == 1602
    skeleton(:,1:cut) = [];
    skeleton(:,end-cut-42:end) = [];
    colony_crop(:,1:cut) = [];
    colony_crop(:,end-cut-42:end) = [];
end

imshow(skeleton<0.5)
circle(x0-cut,y0,rcsr);
circle(x0-cut,y0,rmax);
hold on
plot(x0-cut,y0,'ro')
legend(["rcsr","rmax", "centre"],'FontSize',24)
% exportgraphics(gcf,'skeleton_v2.pdf','ContentType','vector')

%% plot of colony crop
pause(0.1)
figure
imshow(colony_crop<0.5)
circle(x0-cut,y0,rcsr);
circle(x0-cut,y0,rmax);
hold on
plot(x0-cut,y0,'ro')
legend(["rcsr","rmax", "centre"],'FontSize',24)
% exportgraphics(gcf,'colony_crop_v2.pdf','ContentType','vector')

%% plot of original colony fig (a)
pause(0.1)

figure
if size(colony,2) == 1602
    colony(:,1:cut) = [];
    colony(:,end-cut-42:end) = [];
end
imshow(colony<0.5)
% exportgraphics(gcf,'Fig6_binary_50um_s5.pdf','ContentType','vector')
