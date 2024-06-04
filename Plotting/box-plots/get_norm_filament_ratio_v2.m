% this function calculates the filament/colony area ratio of a binary yeast 
% colony image I. 
%
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%
% OUTPUT:
%       
%       filament_norm_ratio: normalised ratio between filamentous area with
%       respect to the whole colony area. 

function filament_norm_ratio = get_norm_filament_ratio_v2(I) 
    
    [x0,y0,~,rcsr] = get_csr_radius_v2(I,3.6);

    % crop colony based on minimum radius 
    M = zeros(size(I,1),size(I,2));
    M(round(y0),round(x0)) = 1;
    R = bwdist(M);
    T = R >= rcsr;    
    colony_outer = T.*I;
    
    % get filamentous to colony area ratio 
    filament_ratio = sum(colony_outer(:))/sum(I(:));
    
    minF = 0.07;  % minimum values of filament_ratio (0.125)
    maxF = 0.18; % maximum values of filament_ratio (0.633)
    
    % normalise to unit interval
    filament_norm_ratio = (filament_ratio-minF)/(maxF-minF);
end