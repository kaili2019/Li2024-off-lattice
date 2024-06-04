% this function calculates the rmin/rmax ratio of a binary yeast colony
% image I.
% 
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%
% OUTPUT:
%       
%       norm_radius_ratio: normalised radius ratio between CSR radius and
%       maximum radius of colony image I. 


function norm_radius_ratio = get_norm_radius_ratio_v2(I)  
    
    [~,~,rmax,rcsr] = get_csr_radius_v2(I,3.6);

    % NOTE: used to be rmin/rmax, but that is increasing trend in the
    % increase in nutrient which is opposite to the other three metrics.
    % Hence, inversed. 
    radius_ratio = rmax/rcsr; 
    
    % min and max radius_ratio to be normalise to unit interval 
    minR = 1.22;
    maxR = 1.45;
    
    % normalise to radius_ratio unit interval
    norm_radius_ratio = (radius_ratio-minR)/(maxR-minR);
end