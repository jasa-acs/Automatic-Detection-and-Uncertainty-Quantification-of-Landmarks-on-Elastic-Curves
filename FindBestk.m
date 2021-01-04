function [post_samp,dsk] = FindBestk(X,maxk,a,b,alpha,numSteps1,numSteps2,var)
%% Produces plot of d_k^2 versus k to pick "optimal" number of lmks k
%% Inputs
% Required:
% X = sample of two-dimensional curves (all curves should be re-sampled to
% the same number of discretization points if they vary)
% maxk = maximum number of landmarks k to consider

% Optional:
% a,b = prior hyperparameters for kappa ~ Gamma(a,b) (default a=1,b=0.01)
% alpha = prior hyperparameter for s ~ Dir(alpha 1) (default alpha=1)
% numSteps1 = # of steps to run MCMC in parallel for all choices of k 
% (default numSteps1=1e6)
% numSteps2 = # of steps to run MCMC for SELECTED value of k (default
% numSteps2=1e6)
% var = variance of proposal density (default var=0.02)
%% Outputs
% post_samp = posterior samples of landmarks theta for user-selected k
% dsk = d_k^2 = average linear reconstruction error for posterior samples
% with k landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check to see if curves are closed
[~,N,~] = size(X);
cl = norm(X(:,N,1)-X(:,1,1)); 
if (cl > 1e-9)
    mink = 1;
else
    mink = 3;
end

% Run MCMC on fixed k model for all possible values of k
parfor i=mink:maxk
    [~,dsk(i)] = ALDfixed(X,i,0,0,a,b,alpha,numSteps1,var)
end

if (cl < 1e-9)
    dsk(1:(mink-1)) = [];
end

% Show elbow plot and have user select optimal value of k
figure(100)
plot(mink:maxk,dsk,'LineWidth',3)
xlabel('k')
ylabel('d_k^2')
set(gca,'FontSize',15)

prompt = 'Which k would you like to select? ';
selectedk = input(prompt);

% Run MCMC again on selected k and produce posterior summaries and plots
[post_samp,~] = ALDfixed(X,selectedk,1,1,a,b,alpha,numSteps2,var);