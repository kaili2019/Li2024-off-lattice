% this function calculate the distance between a simulation
% image and compared to experimental result
%
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%       f_exp: filament area metric of experiment 
%       r_exp: radii ratio metric of experiment
%       s_exp: branch count metric of experiment 
%
% OUTPUT:
%       
%       sim_ss: the simulation summary statistic vector 
%       dist: the distance between the experiment and simulation

function [sim_ss,dist] = calc_summary_statistics_v2(I,f_exp,r_exp,s_exp)
    
    % summary statistics of simulation
    f_I = get_norm_filament_ratio_v2(I);
    r_I = get_norm_radius_ratio_v2(I);
    s_I = get_norm_subbranch_v2(I);
    
    % distance between summary statistics 
    sim_ss = [f_I,r_I,s_I];
    exp_ss = [f_exp,r_exp,s_exp];
    
    % final distance 
    dist = sqrt(sum( ((exp_ss - sim_ss) ./ (1+exp_ss)).^2 ) )/length(sim_ss);
end