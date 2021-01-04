function dst = OpenIntSqDist(X,closest,k,fig)
%% Constructs linear reconstruction thru closest points to landmarks for
%% open curves and computes linear reconstruction error
%% Inputs
% X = two-dimensional curve
% closest = indices of closest discretization points to landmarks
% k = # of landmarks
% fig = display X and linear reconstruction
%% Outputs
% dst = linear reconstruction error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,N] = size(X);

% Construct linear reconstruction
pts = X(:,closest);
intcurve = X(:,1);
for j=1:(length(pts)-1)
    tmpcurve = ReSampleCurve(pts(:,j:j+1),closest(j+1)-closest(j)+1);
    intcurve = [intcurve,tmpcurve(:,2:end)];
end

% Compute linear reconstruction error
qi = curve_to_q(intcurve);
q = curve_to_q(X);
dst = (norm(q(1,:)-qi(1,:))^2)+(norm(q(2,:)-qi(2,:))^2);

% Display original curve and linear reconstruction
if fig==1
    hold on
    axis off
    plot(X(1,:),X(2,:),'LineWidth',3)
    plot(intcurve(1,:),intcurve(2,:),'g','LineWidth',3)
    scatter(pts(1,2:(end-1)),pts(2,2:(end-1)),100,'Filled')
end