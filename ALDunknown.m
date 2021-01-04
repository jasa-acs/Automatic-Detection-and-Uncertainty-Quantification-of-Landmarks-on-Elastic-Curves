function post_samp = ALDunknown(X,lambda,fig,postsumm,a,b,alpha,numSteps,var)
%% Model-based inference of landmark locations and their number 
%% Inputs
% Required:
% X = sample of two-dimensional curves (all curves should be re-sampled to
% the same number of discretization points if they vary)
% lambda = prior hyperparameter for k ~ Shifted Poisson(lambda)
% fig = 1 if want posterior density plots, samples superimposed on curves,
% trace plots, and posterior summaries plotted on curves
% postsumm = 1 if want numerical posterior summaries conditional on chosen
% k (mean, med, MAP, sd, 95% credible intervals)

% Optional:
% a,b = prior hyperparameters for kappa ~ Gamma(a,b) (default a=1,b=0.01)
% alpha = prior hyperparameter for s ~ Dir(alpha 1) (default alpha=1)
% numSteps = # of steps to run RJMCMC (default numSteps=1e5)
% var = variance of proposal density (default var=0.02)

% To use:
% ALDunknown(X,lambda,1,1) - use defaults for all optional parameters
% ALDunknown(X,lambda,1,1,[],0.05,[],1e6,0.03) - use defaults for []
%% Output
% post_samp = posterior samples of landmarks theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup/pre-processing
% d = dimension of curves (=2), N = number of discretization points, 
% M = number of curves in sample
[d,N,M] = size(X);
t = linspace(0,1,N);    % arc-length parameter values (equally spaced)

% Check to see if curves are closed (curves in sample should be all open or
% all closed)
cl = norm(X(:,N,1)-X(:,1,1));   % should be approx 0 if closed

% Re-scale curves to unit Frobenius norm
for m=1:M
    tmp = curve_to_q(X(:,:,m));
    tmp = tmp/sqrt(sum(sum(tmp.^2)));
    X(:,:,m) = q_to_curve(tmp);
end

% Convert to SRVF (after finding reference point for closed curves)
if (cl > 1e-9)
    %% Open curves
    % Compute SRVF
    for m=1:M
        q(:,:,m) = curve_to_q(X(:,:,m));
    end
else
    %% Closed curves
    % Find reference point (pre-processing)
    % Compute first, second derivatives for 1st curve (artificial points
    % inserted before and after initial point of X just for computation)
    [fder,~] = gradient([X(:,end-1,1),X(:,:,1),X(:,2,1)],1/(N+1));
    [sder,~] = gradient(fder,1/(N+1));
    fder(:,[1,end]) = [];   % remove artificially inserted points
    sder(:,[1,end]) = [];
    
    % Compute curvature for first curve
    num = fder(1,:).*sder(2,:)-fder(2,:).*sder(1,:);
    denom = (fder(1,:).^2+fder(2,:).^2).^1.5;
    curv = num./denom;
    
    % Identify point of absolute maximum curvature on 1st curve
    [~,I] = max(abs(curv(1,:)));
    
    % Shift coordinates of all curves to align with reference point
    Xm = X;
    if I~=1     % reference point not equal to initial point
        if M == 1
            Xm(:,end) = [];     % remove repetitive current start/end point
            Xm = ShiftF(Xm,I-1);    
            Xm(:,end+1) = Xm(:,1);  % re-insert to get closed curve
        else
            Xm(:,end,:) = [];
            Xm(:,:,1) = ShiftF(Xm(:,:,1),I-1);  % shift first curve
            q(:,:,1) = curve_to_q(Xm(:,:,1));
            idx_s(1) = I;
            
            % Find shift of remainder of curves in sample which
            % matches first curve optimally (wrt distance on C)
            for m=2:M
                min_d = 5000;   % initialize minimal distance
                for i=1:N
                    tmp = ShiftF(Xm(:,:,m),i-1);
                    tmp_q = curve_to_q(tmp);
                    dst = acos(InnerProd_Q(q(:,:,1),tmp_q));
                    if dst < min_d
                        min_d = dst;
                        idx_s(m) = i;
                    end
                end
                Xm(:,:,m) = ShiftF(Xm(:,:,m),idx_s(m)-1);
                q(:,:,m) = curve_to_q(Xm(:,:,m));   % compute SRVF
            end
            Xm(:,end+1,:) = Xm(:,1,:);  % re-insert to get closed curve
            q(:,end+1,:) = q(:,1,:);
        end
    X = Xm;    
    end
end

%% Generating posterior samples using MCMC
% Set defaults for optional arguments
if ~exist('a','var') || isempty(a), a = 1; end
if ~exist('b','var') || isempty(b), b = 0.01; end
if ~exist('alpha','var') || isempty(alpha), alpha = 1; end
if ~exist('numSteps','var') || isempty(numSteps), numSteps = 1e5; end
if ~exist('var','var') || isempty(var), var = 0.02; end

% Other MCMC settings
burnIn = 0.1*numSteps;      % number of burn-in samples (adjust if slow-mixing)
thin = 100;                 % thinning number
AP(1) = 1;                  % initialize acceptance probability
accept(1) = 1;              % initialize acceptance indicator

% Establish minimum number of landmarks and initialize k
if (cl > 1e-9)
    k_min = 1;
else
    k_min = 3;
end

k(1) = poissrnd(lambda)+k_min;  % initialize k

% Initialize birth/death/stay probabilities 
p(:,1) = MoveProb(k(1),k_min);

% Size of set s corresponding to theta
if (cl > 1e-9)
    s_dim(1) = k(1)+1;
else
    s_dim(1) = k(1);
end

% Concentration parameter for s (alpha)
conc{1} = alpha*ones(1,s_dim(1));

% Initialize theta/s
theta{1} = ((1:k(1))/(k(1)+1))';
if (cl > 1e-9)
    S{1} = diff([0;theta{1};1]);    
else
    S{1} = [diff(theta{1});mod(theta{1}(1)-theta{1}(end),1)];
end

% Identify closest point in discretization to theta values
[~,closest] = min(abs(repmat(t,k(1),1)-repmat(theta{1},1,N)),[],2);
if (cl > 1e-9)
    closest = [1;closest;N];    % include start/end for open curves
end

% Compute cumulative linear reconstruction error
dst(1) = 0;
if (cl > 1e-9)
    for m=1:M
        dst(1) = dst(1) + OpenIntSqDist(X(:,:,m),closest,k(1),0);
    end
else
    for m=1:M
        dst(1) = dst(1) + ClosedIntSqDist(X(:,:,m),closest,k(1),0);
    end
end

% Log prior, log likelihood, log posterior
lprior(1) = log(gamma(sum(conc{1})))-sum(log(gamma(conc{1})))+sum((conc{1}-1)'.*log(S{1}))-lambda+(k(1)-k_min)*log(lambda)-log(factorial(k(1)-k_min));
llik(1) = -N*M*log(pi)+log(gamma(a+N*M))+a*log(b)-log(gamma(a))-(a+N*M)*log(b+dst(1));
lpost(1) = lprior(1)+llik(1);

% Run chain
for n=2:numSteps
    clear propS
    
    % Determine move type
    j = rand;
    if j <= p(1,n-1)    % birth
        prop_k = k(n-1)+1;
        mov(n) = 1;
    elseif p(1,n-1) < j && j <= (p(1,n-1)+p(2,n-1))    % death
        prop_k = k(n-1)-1;
        mov(n) = 2;
    else                % stay
        prop_k = k(n-1);
        mov(n) = 3;
    end
    
    % Updated move probabilities
    prop_p = MoveProb(prop_k,k_min);
    
    % Generate proposal vector theta^*
    prop = theta{n-1};
    
    if mov(n) == 1          % birth
        u = rand;
        prop = sort([prop;u]);
    elseif mov(n) == 2      % death
        idx_d = randi([1,k(n-1)]);
        prop(idx_d) = [];
    else                    % stay
        idx_s = randi([1,k(n-1)]);
        tmp = prop(idx_s);
        update = normrnd(tmp,var);
        if (cl > 1e-9)
            prop(idx_s) = update;
        else
            prop(idx_s) = mod(update,1);
            
            % Account for proposals below 0 or above 1 to re-sort proposal
            if update < 0
                prop = (ShiftF(prop',1))';
            elseif update > 1
                prop = (ShiftF(prop',prop_k-1))';
            end
        end
    end
        
    % Compute corresponding set s from theta^*
    if (cl > 1e-9)
        propS = diff([0;prop;1]);
    else
        propS = [diff(prop);mod(prop(1)-prop(end),1)];
    end
    
    % Find closest discretization points to proposed landmarks
    [~,closest] = min(abs(repmat(t,prop_k,1)-repmat(prop,1,N)),[],2);
    if (cl > 1e-9)
        closest = [1;closest;N];
    end
    
    % Compute linear reconstruction error and log
    % prior/likelihood/posterior for proposal, and accept/reject via M-H
    propdst = 0;
    
    % Automatically reject non-unique/invalid proposed landmark values
    if ((cl > 1e-9) && ( length(unique(closest))<(prop_k+2) || all(propS>0)==0 || any(prop)<0 || any(prop)>1 ))
        AP(n) = 0;
        accept(n) = 0;
        theta{n} = theta{n-1};
        S{n} = S{n-1};
        lprior(n) = lprior(n-1);
        llik(n) = llik(n-1);
        lpost(n) = lpost(n-1);
        dst(n) = dst(n-1);   
        s_dim(n) = s_dim(n-1);
        conc{n} = conc{n-1};
        k(n) = k(n-1);
        p(:,n) = p(:,n-1);
    elseif ((cl < 1e-9) && ( length(unique(closest))<prop_k || closest(1)==closest(prop_k)-(N-1) || all(propS>0)==0 ))
        AP(n) = 0;
        accept(n) = 0;
        theta{n} = theta{n-1};
        S{n} = S{n-1};
        lprior(n) = lprior(n-1);
        llik(n) = llik(n-1);
        lpost(n) = lpost(n-1);
        dst(n) = dst(n-1);
        s_dim(n) = s_dim(n-1);
        conc{n} = conc{n-1};
        k(n) = k(n-1);
        p(:,n) = p(:,n-1);
    else
        % Compute proposal linear reconstruction error
        if (cl > 1e-9)
            for m=1:M
                propdst = propdst + OpenIntSqDist(X(:,:,m),closest,prop_k,0);
            end
        else
            for m=1:M
                propdst = propdst + ClosedIntSqDist(X(:,:,m),closest,prop_k,0);
            end
        end
        
        % Size of set s corresponding to theta
        if (cl > 1e-9)
            s_dim_prop = prop_k+1;
        else
            s_dim_prop = prop_k;
        end
        
        % Concentration parameter for s (alpha)
        propconc = alpha*ones(1,s_dim_prop);
        
        % Proposal log prior/likelihood/posterior
        proplprior = log(gamma(sum(propconc)))-sum(log(gamma(propconc)))+sum((propconc-1)'.*log(propS))-lambda+(prop_k-k_min)*log(lambda)-log(factorial(prop_k-k_min));
        propllik = -N*M*log(pi)+log(gamma(a+N*M))+a*log(b)-log(gamma(a))-(a+N*M)*log(b+propdst);
        proplpost = proplprior+propllik;
        
        % Metropolis ratio and acceptance probability
        if mov(n) == 1      % birth move
            MH(n) = exp(proplpost+log(prop_p(2))-lpost(n-1)-log(prop_k)-log(p(1,n-1)));
        elseif mov(n) == 2  % death move
            MH(n) = exp(proplpost+log(prop_p(1))+log(k(n-1))-lpost(n-1)-log(p(2,n-1)));
        else
            MH(n) = exp(proplpost-lpost(n-1));
        end

        AP(n) = min([1,MH(n)]);
        
        if rand < AP(n)     % accept proposal
            accept(n) = 1;
            theta{n} = prop;
            S{n} = propS;
            lprior(n) = proplprior;
            llik(n) = propllik;
            lpost(n) = proplpost;
            dst(n) = propdst;
            s_dim(n) = s_dim_prop;
            conc{n} = propconc;
            k(n) = prop_k;
            p(:,n) = prop_p;
        else                % reject proposal
            accept(n) = 0;
            theta{n} = theta{n-1};
            S{n} = S{n-1};
            lprior(n) = lprior(n-1);
            llik(n) = llik(n-1);
            lpost(n) = lpost(n-1);
            dst(n) = dst(n-1);
            s_dim(n) = s_dim(n-1);
            conc{n} = conc{n-1};
            k(n) = k(n-1);
            p(:,n) = p(:,n-1);
        end
    end
    
    % Print every 1000th step
    if mod(n,1000)==0
        n
    end
end

acc_rate = sum(accept)/numSteps;            % acceptance rate
post_samp = theta(:,(burnIn+1):thin:end);   % remove burn-in and thin

% Group posterior theta samples of the same dimension together
post_size = cellfun(@length,post_samp);     % records dim of post. samples
rangek = max(post_size)-min(post_size);

for i=k_min:(max(post_size))
    theta_k{i} = post_samp(post_size==i);
end

% Resample original curves to a much larger size for visualizing posteriors
N1 = 10000;
for m=1:M
    for i=1:d
        X1(i,:,m) = interp1(t,X(i,:,m),linspace(0,1,N1));
    end
end

%% Histogram for posterior distribution of k
if (fig==1)
    histogram(post_size,'Normalization','probability')
    xlabel('k')
    ylabel('Probability')
    set(gca,'FontSize',20)
end

%% Posterior for theta given a chosen value in posterior of k
post_k = min(post_size):max(post_size)
prompt = 'Which k in the posterior would you like to select? (must be in previous list) ';
k_fig = input(prompt);

theta_k_tmp = cell2mat(theta_k{k_fig});

% Find landmark locations on corresponding curves
for i=1:length(theta_k_tmp)
    for m=1:M
        [~,closest] = min(abs(repmat(linspace(0,1,N1),k_fig,1)-repmat(theta_k_tmp(:,i),1,N1)),[],2);
        land(:,:,i,m) = X1(:,closest,m);
    end
end

% Post-processing of posterior for closed curves
if (cl < 1e-9)
    tmpsamp = theta_k_tmp;
    for i=2:length(tmpsamp)
        for j=1:k_fig
            shift = ShiftF(tmpsamp(:,i)',j)';
            for l=1:k_fig
                tmp = theta_k_tmp(l,1)-shift(l);
                d_c(l) = min([abs(tmp),abs(tmp+1),abs(tmp-1)],[],2);
            end
            d(j) = sum(d_c);
        end
        [~,ind] = min(d);
        for m=1:M
            land(:,:,i,m)=ShiftF(land(:,:,i,m),ind);
        end
        theta_k_tmp(:,i) = ShiftF(tmpsamp(:,i)',ind)';
    end
end

%% Compute posterior summaries conditional on chosen k
% Mean/median/sd/credible intervals
if (cl > 1e-9)
    meantheta = mean(theta_k_tmp,2);
    medtheta = median(theta_k_tmp,2);
    sdtheta = std(theta_k_tmp,[],2);
    thetacred = prctile(theta_k_tmp,[2.5,97.5],2);
else
    % For mean, median, sd - need to manually adjust values of first 
    % component of theta to account for circular domain, i.e., if values of 
    % theta straddle 0-1 boundary, "unmod" them by adding 1 to all values 
    % which are close to 0
    
    % Example: samples are theta = [0.92,0.94,0.98,0.03,0.07];
    % newtheta = theta + (theta<0.4) will yield [0.92,0.94,0.98,1.03,1.07]
    % Then take mean/median/sd of newtheta and for mean/median, mod by 1
    adjsamp1 = theta_k_tmp(1,:)+(theta_k_tmp(1,:)<0.4); % adjust 0.4 by case
    meantheta(1) = mod(mean(adjsamp1),1);
    medtheta(1) = mod(median(adjsamp1),1);
    sdtheta(1) = std(adjsamp1);
    thetacred(1,:) = mod(prctile(adjsamp1,[2.5,97.5]),1);

    meantheta(2:k_fig) = mean(theta_k_tmp(2:k_fig,:),2)';
    medtheta(2:k_fig) = median(theta_k_tmp(2:k_fig,:),2)';
    sdtheta(2:k_fig) = std(theta_k_tmp(2:k_fig,:),[],2)';
    thetacred(2:k_fig,:) = mod(prctile(theta_k_tmp(2:k_fig,:),[2.5,97.5],2),1);
    
    meantheta = meantheta';
    medtheta = medtheta';
    sdtheta = sdtheta';
end

if (postsumm == 1)
    varNames = {'Landmark','Mean' 'Median' 'SD' 'CI_95'};
    table((1:k_fig)',meantheta,medtheta,sdtheta,thetacred,'VariableNames',varNames)
end

%% Figures
if (fig == 1)
    % Closest points on curves to mean/median/MAP/credible interval estimates
    [~,closest] = min(abs(repmat(linspace(0,1,N1),k_fig,1)-repmat(meantheta,1,N1)),[],2);
    meanpts = X1(:,closest,:);
    
    [~,closest] = min(abs(repmat(linspace(0,1,N1),k_fig,1)-repmat(medtheta,1,N1)),[],2);
    medpts = X1(:,closest,:);
    
    for i=1:size(thetacred,2)
        [~,closest] = min(abs(repmat(linspace(0,1,N1),k_fig,1)-repmat(thetacred(:,i),1,N1)),[],2);
        ciend(:,:,:,i) = X1(:,closest,:);
    end

    % Marginal posterior densities of theta
    for i=1:k_fig
        figure(i)
        [x1,x2] = ksdensity(theta_k_tmp(i,:),'function','pdf');
        plot(x2,x1,'LineWidth',3)
        str = sprintf('Landmark %i',i);
        title(str);
        set(gca,'FontSize',20)
    end 
    
    % Original curves with posterior landmark samples superimposed
    for m=1:M
        figure(k_fig+m)
        hold on
        plot(X1(1,:,m),X1(2,:,m),'LineWidth',3)
        for i=1:k_fig
            scatter(land(1,i,:,m),land(2,i,:,m),90)
        end
        axis off equal
    end
    
    % Posterior summaries superimposed on original curves
    for m=1:M
        figure(k_fig+M+m+1);
        hold on
        axis off equal
        plot(X1(1,:,m),X1(2,:,m),'LineWidth',3)
        for i=1:size(land,2)
            scatter(ciend(1,i,m,1:2),ciend(2,i,m,1:2),150,'s','filled')
        end
        ax = gca;
        ax.ColorOrderIndex = 2;
        for i=1:size(land,2)
            scatter(meanpts(1,i,m),meanpts(2,i,m),150,'o','LineWidth',2)
        end
        ax = gca;
        ax.ColorOrderIndex = 2;
        for i=1:size(land,2)
            scatter(medpts(1,i,m),medpts(2,i,m),150,'*','LineWidth',2)
        end
    end
end