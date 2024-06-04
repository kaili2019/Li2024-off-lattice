% this function finds the complete spatial randomness (CSR) radius of a
% given binary yeast colony rcsr. It also outputs the center of the colony
% as well as the maximum radius. 
%
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%       del: typically set to 4. See Binder et al. (2015). 
%
% OUTPUT:
%       
%       x0: the x-coordinates of the centre
%       y0: the y-coordinates of the centre
%       rmax: maximum raidus of colony 
%       rcsr: the complete spatial randomness radius  

function [x0,y0,rmax,rcsr] = get_csr_radius_v2(I,del)
    
    % get centre of colony
    measurements = regionprops(I, 'Centroid');
    Image_centre = measurements.Centroid;
    
    % find maximum radius 
    dist = bwdist(I);
    
    [y_max,x_max]=find(dist==0);
    
    xc_max = x_max-Image_centre(1);
    yc_max = y_max-Image_centre(2);
    
    E_dist_max = sqrt(xc_max.^2+yc_max.^2);
    
    rmax = max(E_dist_max);
    
    % calculate mean field density. See Binder et al (2015). 
    
    [Y,X] = find(I);
    x0 = Image_centre(1);
    y0 = Image_centre(2);
    
    rho = sum(I(:))/(pi*rmax^2);
    L = round(rmax/del);
    fr = nan(L,1);
    
    for i = 1:L
        n_mi = sum(sum(((X-x0).^2+(Y-y0).^2 <= (del*(i-1))^2)));
        n_i = sum(sum(((X-x0).^2+(Y-y0).^2 <= (del*i)^2)));
        cr = n_i - n_mi;
        fr(i) = cr/(pi*rho*del^2*(2*i-1));
    end

    [~, min_idx] = min(abs(1-fr));

    % calculate the CSR radius 
    rcsr = min_idx*del;
end
