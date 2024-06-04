% this function take a set of parameter values and simulates the colony to
% the same area as the experiment upto 5% difference. 
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
%       upr_area: upper bound of area within 5% of experimental area
%       lwr_area: lower bound of area within 5% of experimental area
%       exp_area: experimental colony pixel area
%
% OUTPUT:
%       I: a binary colony of the simulation
%       file_name: the file name associated to the simualtion including the
%       information related to the particular simulation


function [I,file_name] = run_off_lattice(N,Telong,p2sProb,s2pProb,pc,pa,upr_area,lwr_area,exp_area)

    % run simualtion 
    [a2,b2,b_el2,pos] = simulate_colony(N,Telong,p2sProb,s2pProb,pc,pa);
    
    % simulates colony with area correction 
    
    % calculate (x,y) coordinates of cells within colony
    [x,y] = calc_colony_coordinates(a2,b2,b_el2,pos,N);
    
    % saving colony as tiff
    file_name = "Nfinal="+N+"_N="+N+"_Telong="+Telong+"_p2sProb="+p2sProb+"_s2pProb"...
                +s2pProb+"_pc="+pc+"_pa="+pa; 
    
    % save colony as matrix
    [I_hist,~,~] = histcounts2(x+700, y+900,0:4:1400,0:4:1800);
    I = imresize(I_hist>0.5,4);
    I = ~bwareaopen(~I, 10); % fill gaps in colony that has 10 or less pixels

    I_area = sum(I(:));
    
    % check area of experiment is below simulations area
    if exp_area < I_area
        area_guess = exp_area/I_area - 0.1;
        same_area = 0;
    else
        same_area = 1;
        Nguess = N;
        disp("Simultaed colony area is smaller than experiment!")
    end
    

    % area correction 
    while ~same_area
        
        Nguess = round(area_guess*N);
    
        % calculate (x,y) coordinates of cells within colony
        [x,y] = calc_colony_coordinates(a2,b2,b_el2,pos,Nguess);
        
        % saving colony as tiff
        file_name = "Nfinal="+N+"_N="+Nguess+"_Telong="+Telong+"_p2sProb="+p2sProb+"_s2pProb"...
                    +s2pProb+"_pc="+pc+"_pa="+pa; 
        
        % save colony as matrix
        [I_hist,~,~] = histcounts2(x+700, y+900,0:4:1400,0:4:1800);
        I = imresize(I_hist>0.5,4);
        I = ~bwareaopen(~I, 10); % fill gaps in colony that has 10 or less pixels
        
        I_area = sum(I(:));
        
        if I_area > lwr_area && I_area < upr_area
            same_area = 1;
        end
        
        if area_guess+0.025 < 1
            area_guess = area_guess+0.025;
        else
            area_guess = 1;
            same_area = 1;
        end
    
    end

end