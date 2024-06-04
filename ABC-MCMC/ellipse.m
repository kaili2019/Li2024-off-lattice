% This function takes in 
% 
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       
%
% OUTPUT:
%       
%       

function [x,y] = ellipse(x0, y0,a,b,theta)
    phi = linspace(0,2*pi-pi/6,12);
    x = x0 + a*cos(phi)*cos(theta) - b*sin(phi)*sin(theta);
    y = y0 + b*sin(phi)*cos(theta) + a*cos(phi)*sin(theta);
end