function [x,y] = calc_colony_coordinates(a2,b2,b_el2,pos,N)

    % calculating the x,y coordinates of cells 
    phi = linspace(0,2*pi-pi/6,12);
    
    % start printing from this cell onwards 
    start = 1;
    stop = N;
    
    x0 = pos(start:stop,1); 
    y0 = pos(start:stop,2);
    theta = pos(start:stop,3);
    a = pos(start:stop,4);
    
    % calculate length b based on length a
    b = (a==a2)*b2 + (a~=a2)*b_el2;
    
    x = x0 + a*cos(phi).*cos(theta) - b*sin(phi).*sin(theta);
    y = y0 + b*sin(phi).*cos(theta) + a*cos(phi).*sin(theta);

end