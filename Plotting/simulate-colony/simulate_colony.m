% this function simulate colonies given a set of parameter values
% 
% AUTHOR: KAI LI
% 
% Date: 22 Mar 2024
% 
% INPUT: 
%       N: run colony upto this amount of cells
%       Telong: Proportion of cells to be sated before switching on 
%               pseudohyphal cells probabilities
%       p2sProb: Probability of pseudohyphal cells transitioning to 
%               sated cells
%       s2pProb: Probability of sated cells transitioning to pseudohyphal 
%               cells
%       pc: Probability of second daughter cell being sated given current 
%               cell is pseudohyphal
%       pa: Probability of choosing a sated cell to proliferate
%
% OUTPUT:
%       a2: width of sated cell
%       b2: height of sated cell
%       b_el3: height of pseudohyphal cell
%       pos: a matrix of coordinates of the colony. column 1 and 2 are the
%       (x,y) coordinates, column 3 is the beta angle relative to the
%       x-axis, and column 4 is the a vector to identify the type of the
%       cell by saving the width of the cell. 

function [a2,b2,b_el2,pos] = simulate_colony(N,Telong,p2sProb,s2pProb,pc,pa)

    L = 323/500; % 1 pixel per um
    
    a = 4.2*L;
    b = 3.0*L;
    a2 = a;
    b2 = b;
    
    a_el = 6.7*L;
    b_el = 1.9*L;
    a_el2 = a_el;
    b_el2 = b_el;
     
    Beta = 7*pi/16;
    betap = pi/16;
    
    % Total number of nutrients available 
    % N = 30000;
    
    % time to transition to pseudohyphal mode 
    % Telong = 0.75;
    % probability of transitioning from pesudohyphal to sated
    % p2sProb = 0.3; 
    % probability of transitioning from sated to pesudohyphal 
    % s2pProb = 0.3;
    % probability of second cell being pseudohyphal 
    % pc = 0.05;
    % probability of choosing pseudhyphal cell
    % pa = 0.05;
    % number to first implementing pa rule
    Npseudo = (Telong+0.01)*N;
    
    % intialise angles 
    theta = 0; % angle relative to x-axis 
    angle = [Beta, -Beta, pi+Beta, pi-Beta]; % angles of proliferation site occurs
    ang_sated = angle;
    % thetap = 0; % angle relative to x-axis for enlongated cell
    ang_elong = [betap, -betap, betap, -betap];
    
    center = 0; % center of first cell
    
    % number of sides on polygon for approximating ellipse 
    len = 12;
    
    % intialise storage of x,y-coordinate of cells edge
    ex = NaN(N*len+N-1,1);
    ey = NaN(N*len+N-1,1);
    
    % store x,y-coordinate of initial cell edges
    [ex(1:len,:),ey(1:len,:)] = ellipse(center, center,a,b,0);
    
    % declare functions
    
    % angle of proliferation wrt the prolieration cell function
    phib = @(beta,ar,br) beta + (atan((ar-br)*tan(beta)))/(br+ar*tan(beta)^2);
    
    % x coordinate on border of ellipse function
    x = @(x0,phi,theta,ar,br) x0 + ar*cos(phi)*cos(theta) ...
        - br*sin(phi)*sin(theta);
                    
    % y coordinate on border of ellipse function
    y = @(y0,phi,theta,ar,br) y0 + ar*cos(phi)*sin(theta) ...
        + br*sin(phi)*cos(theta);
    
    % store x,y-coordinates of centre of cell and  proliferation site angle
    pos = nan(N,4);
    
    % store initial cell
    pos(1,1) = center;
    pos(1,2) = center;
    pos(1,3) = theta;
    pos(1,4) = a;
    
    % time to introduce elongated cells
    t = round(Telong*N);
    
    keep = NaN(N,1);
    keep_pseudo = zeros(N,1);
    
    % initialise the initial cell to have zero prolierations done
    keep(1) = 1;
    keep_pseudo(1) = 0;

    tic
    for i = 2:N % iterate over all cells
    
        % select random cell for proliferation that CAN proliferate
        id_array = find(keep(1:i) < 4);
        id = id_array(randi(length(id_array)));
    
        % introduce elongated cell
        if i > t 
            prob = rand; % choose an random probability 
            prob2 = rand;
            prob3 = rand;
    
            if prob3 >= pa && i >= Npseudo
                id_array = find(keep(1:i) < 4 & pos(1:i,4) == a_el2);
                id = id_array(randi(length(id_array)));
            elseif prob3 < pa && i < Npseudo
                id_array = find(keep(1:i) < 4 & pos(1:i,4) == a2);
                id = id_array(randi(length(id_array)));
            end
    
            if prob2 > pc && keep_pseudo(id) ~= 0 && prob > p2sProb 
                id_array = find(keep(1:i) < 4 & keep_pseudo(1:i) == 0);
                id = id_array(randi(length(id_array)));
            end
    
            % sated to pseudohyphal 
            if pos(id,4) == a2 && prob < s2pProb
                a = a_el;
                b = b_el;
                angle = ang_elong;
            % sated to sated 
            elseif pos(id,4) == a2 && prob > s2pProb
                a = a2;
                b = b2;
                angle = ang_sated;
            % pseudohyphal to pseudohyphal (must first produce a pseudohyphal
            % cell)
            elseif pos(id,4) == a_el2 && keep_pseudo(id) == 0 
                a = a_el;
                b = b_el;
                angle = ang_elong;
            % pseudohyphal to sated
            elseif pos(id,4) == a_el2 && keep_pseudo(id) ~= 0 && prob < p2sProb 
                a = a2;
                b = b2;
                angle = ang_sated;
            % pseudohyphal to second pseudohyphal cell
            elseif pos(id,4) == a_el2 && keep_pseudo(id) ~= 0 && prob > p2sProb 
                a = a_el;
                b = b_el;
                angle = ang_elong;
            else
                disp("Erorr caught found in if statement")
                break
            end
        end
    
        theta = pos(id,3); % get angle of current cell
        x0 = pos(id,1); % get x-coordinate of current cell
        y0 = pos(id,2); % get y-coordinate of current cell
        
        Beta = angle(randi(length(angle))); % get angle of new proliferation site
        
        % check if to introduce elongated cell
        if id <= t
            X = x(x0,phib(Beta,a2,b2),theta,a2,b2) + a*cos(Beta + theta);
            Y = y(y0,phib(Beta,a2,b2),theta,a2,b2) + a*sin(Beta + theta);
        else
            if pos(id,4) == a2 % if sated 
                X = x(x0,phib(Beta,a2,b2),theta,a2,b2) + a*cos(Beta + theta);
                Y = y(y0,phib(Beta,a2,b2),theta,a2,b2) + a*sin(Beta + theta);
            else % if pseudohyphal 
                X = x(x0,phib(Beta,a_el,b_el),theta,a_el,b_el) + a*cos(Beta + theta);
                Y = y(y0,phib(Beta,a_el,b_el),theta,a_el,b_el) + a*sin(Beta + theta);
            end
        end
        
        n_id = find(sqrt((X-pos(1:i,1)).^2+(Y-pos(1:i,2)).^2) < 3*a);
        
        % create start and end index for x,y-coordinate of edge of polygon
        % ellipse
        id_start = (n_id-1)*len+n_id;
        
        n_id_array =  id_start' + (0:len)';
        n_id_array = n_id_array(:);
            
        % check if center is in cell
        while inpolygon(X,Y,ex(n_id_array),ey(n_id_array)) 
            % select random cell for proliferation that CAN proliferate
            id_array = find(keep(1:i) < 4);
            
            % get id of new prolifferation cell
            id = id_array(randi(length(id_array)));
    
           % introduce elongated cell
            if i > t 
                prob = rand; % choose an random probability 
                prob2 = rand;
                prob3 = rand;
        
                if prob3 >= pa && i >= Npseudo
                    id_array = find(keep(1:i) < 4 & pos(1:i,4) == a_el2);
                    id = id_array(randi(length(id_array)));
                elseif prob3 < pa && i < Npseudo
                    id_array = find(keep(1:i) < 4 & pos(1:i,4) == a2);
                    id = id_array(randi(length(id_array)));
                end
    
                if prob2 > pc && keep_pseudo(id) ~= 0 && prob > p2sProb
                    id_array = find(keep(1:i) < 4 & keep_pseudo(1:i) == 0);
                    id = id_array(randi(length(id_array)));
                end
    
                % sated to pseudohyphal 
                if pos(id,4) == a2 && prob < s2pProb
                    a = a_el;
                    b = b_el;
                    angle = ang_elong;
                    % sated to sated 
                elseif pos(id,4) == a2 && prob > s2pProb
                    a = a2;
                    b = b2;
                    angle = ang_sated;
                    % pseudohyphal to pseudohyphal (must first produce a pseudohyphal
                    % cell)
                elseif pos(id,4) == a_el2 && keep_pseudo(id) == 0 
                    a = a_el;
                    b = b_el;
                    angle = ang_elong;
                    % pseudohyphal to sated
                elseif pos(id,4) == a_el2 && keep_pseudo(id) ~= 0 && prob < p2sProb 
                    a = a2;
                    b = b2;
                    angle = ang_sated;
                    % pseudohyphal to second pseudohyphal cell
                elseif pos(id,4) == a_el2 && keep_pseudo(id) ~= 0 && prob > p2sProb 
                    a = a_el;
                    b = b_el;
                    angle = ang_elong;
                else
                    disp("Erorr caught found in if statement")
                    break
                end
            end
    
            
            theta = pos(id,3);
            Beta = angle(randi(length(angle)));
            x0 = pos(id,1);
            y0 = pos(id,2);
            
            % check if to introduce elongated cell
            if id <= t
                X = x(x0,phib(Beta,a2,b2),theta,a2,b2) + a*cos(Beta + theta);
                Y = y(y0,phib(Beta,a2,b2),theta,a2,b2) + a*sin(Beta + theta);
            else
                if pos(id,4) == a2 
                    X = x(x0,phib(Beta,a2,b2),theta,a2,b2) + a*cos(Beta + theta);
                    Y = y(y0,phib(Beta,a2,b2),theta,a2,b2) + a*sin(Beta + theta);
                else
                    X = x(x0,phib(Beta,a_el,b_el),theta,a_el,b_el) + a*cos(Beta + theta);
                    Y = y(y0,phib(Beta,a_el,b_el),theta,a_el,b_el) + a*sin(Beta + theta);
                end
            end
            
            n_id = find(sqrt((X-pos(1:i,1)).^2+(Y-pos(1:i,2)).^2) < 3*a);
            
            id_start = (n_id-1)*len+n_id;
    
            n_id_array =  id_start' + (0:len)';
            n_id_array = n_id_array(:);
            
            if i < t
                keep(id) = keep(id)+1;
            end
        end
        
        if i >= t
            keep(id) = keep(id)+1;
        end
        
        % update to new proliferation cell 
        [ex((i-1)*len+i:i*len+(i-1),1),ey((i-1)*len+i:i*len+(i-1),1)] = ...
            ellipse(X, Y,a,b,theta+Beta);
        
        % store new cell info for later use 
        pos(i,1) = X; 
        pos(i,2) = Y;
        pos(i,3) = Beta + theta;
        pos(i,4) = a;
        
        keep(i) = 1;
        
        if a == a_el2 
            keep_pseudo(id) = keep_pseudo(id)+1;
        end
        
%         if mod(i,10000) == 0
%             fprintf("%0.2f precentage simulated took %0.2f seconds" + ...
%                 " with N = %d\n", i/N*100, toc, i);
%             time_vec(count,1) = toc;
%             time_vec(count,2) = i;
%             count = count + 1;
%         end
    end  

end