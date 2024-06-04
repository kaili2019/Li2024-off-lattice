% this function calculates the number of subbranches of a colony based on
% the skeletonized image. 
% 
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%
% OUTPUT:
%       
%       norm_subbranch: normalised subbranch count from the skeletonization
%       method. 

function norm_subbranch = get_norm_subbranch_v2(I)

    % get csr radius
%     [x0,y0,rmin,~] = get_radii(I);
    
    [x0,y0,~,rcsr] = get_csr_radius_v2(I,4);

    % get a matrix to crop out the centre of the colony based on csr radius
    M = zeros(size(I,1),size(I,2));
    M(round(y0),round(x0)) = 1;
    R = bwdist(M);
    C = R >= rcsr;
    
    % crop a quarter of the colony
    T = zeros(size(I,1),size(I,2));
    T(1:round(y0),1:round(x0)) = 1;  

    % skeletonize
    skeleton = bwmorph(I.*T.*C, 'skel', Inf);

    % calculate sub branch number 
    It = logical(skeleton);
    [i,~] = find(bwmorph(It,'endpoints'));
    
    % This is normalising assuming the least number of branch count is minS 
    % while the most is maxS. From the experimental data the largest subbranch 
    % count recorded was 597 and minimum is 54.
    minS = 25;
    maxS = 600;
        
    subbranch = length(i);

    % normalise subbranch count unit interval
    norm_subbranch = (subbranch-minS)/(maxS-minS);
end