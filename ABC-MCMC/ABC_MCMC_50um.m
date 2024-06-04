% This script implement the ABC-MCMC on inferring the five parameter values 
% of the off-lattice model. In particular this infers the parameters values
% of the AWRI 50um.
% 
% Kai Li
% 18 April 2024

clear 
clc
close all

rng(20) % set random seed for HPC use 

load("cov_50_v15.mat")

disp(cov_50_v15)

% rescale covariance matrix by 1/3 to reduce step size
cov_50_v15 = cov_50_v15/3;

disp(cov_50_v15)

cc = 1; % count number of accepted 

% beta distribution parameters 
A = [5, 2, 2, 2, 2]; 
B = [2, 2, 2, 5, 5];

% MC1. Initialize theta_1
Nsims = 800; % number of iterations
theta = nan(Nsims,5); % matrix of parameter values 
theta(1,:) = [0.65,0.5,0.3,0.08,0.1] + 0.01*rand; % intial guess 

% time to transition to pseudohyphal mode 
Telong = theta(1,1);
% probability of transitioning from pesudohyphal to sated
p2sProb = theta(1,2); 
% probability of transitioning from sated to pesudohyphal 
s2pProb = theta(1,3); 
% probability of second cell being pseudohyphal 
pc = theta(1,4); 
% probability of choosing pseudhyphal cell
pa = theta(1,5); 

% threshold 
eps = 0.03;

%% experimental image summary statistic

load("ss_mean_50.mat")

f_exp = ss_mean_50(1);
r_exp = ss_mean_50(2);
s_exp = ss_mean_50(3);

exp_area =  573920; % from the average of 50um colonies

% run simulations upto 5% of the experimental colony size
lwr_area = exp_area*0.95;
upr_area = exp_area*1.05;

% initialise array to store information from ABC-MCMC
samples = cell(Nsims,6);
Edist = nan(Nsims,1); 

% total number of nutrients available 
N = 45000;

%%
tStart = tic;
for i = 2:Nsims
    tic
    % MC2. Generate a candidate value theta^* ~ q(theta|theta_i) where q(.) 
    % assumed to be multinormal distributed. 
    thetaStar = mvnrnd(theta(i-1,:),cov_50_v15);

    % MC3. Generate a data set x^* ~ f(x|thetaStar) 
    Telong = thetaStar(1);
    p2sProb = thetaStar(2);
    s2pProb = thetaStar(3); 
    pc = thetaStar(4); 
    pa = thetaStar(5); 
    
    % check if outside bound, if yes, then reject
    if sum((thetaStar >= 0) & (thetaStar <= 1)) == 5 
        % run simulation 
        [I,file_name] = run_off_lattice(N,Telong,p2sProb,s2pProb,pc,pa,upr_area,lwr_area,exp_area);
        
        [sim_ss,Edist(i)] = calc_summary_statistics_v2(I,f_exp,r_exp,s_exp);
    
        % MC4. Set theta_1 = thetaStar with probability alpha
        numerator = prod(betapdf(thetaStar,A,B));
        denominator = prod(betapdf(theta(i-1,:),A,B));

        % accept samples with probability alpha
        alpha = min([1,(numerator/denominator)*(Edist(i) <= eps)]);
    else
        alpha = 0;
    end

    if rand < alpha % accepts
        theta(i,:) = thetaStar;
        fprintf("SAMPLE ACCEPTED. %d samples have been accepted with alpha = %0.2f.\n",cc,alpha);
        cc = cc+1;
    else % reject 
        theta(i,:) = theta(i-1,:);
    end
    
    % save information from chain 
    samples{i,1} = I;
    samples{i,2} = file_name;
    samples{i,3} = sim_ss;
    samples{i,4} = Edist(i);
    samples{i,5} = [Telong,p2sProb,s2pProb,pc,pa];
    samples{i,6} = theta(i,:);

    fprintf("%d simulated completed and took %0.2f seconds with %d simulations left.\n", i-1, toc, Nsims-i);
    fprintf("The distance was %0.3f\n",Edist(i))
    fprintf("Telong = %0.3f, p2sProb = %0.3f, s2pProb = %0.3f, pc = %0.3f, pa = %0.3f\n\n",...
            Telong,p2sProb,s2pProb,pc,pa);

    % note: used MATLAB 2020 version due to HPC purposes

    if mod(i,100) == 0 
        save("simulations/theta v1 i="+i+" "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","theta")
        save("simulations/dist v1 i="+i+" "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","Edist")
        save("simulations/samples v1 i="+i+" "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
    end
    
end

tEnd = toc(tStart);

fprintf("total time taken is %0.2f seconds\n",tEnd);

%%

save("final_theta "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","theta")
save("final_dist "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","Edist")
save("final_samples "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
