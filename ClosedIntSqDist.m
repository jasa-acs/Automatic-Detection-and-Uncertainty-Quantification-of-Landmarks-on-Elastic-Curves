function dst = ClosedIntSqDist(X,closest,k,fig)
%% Constructs linear reconstruction thru closest points to landmarks for
%% closed curves and computes linear reconstruction error
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
tmpcurve = ReSampleCurve(pts(:,1:2),closest(2)-closest(1)+1);
intcurve = tmpcurve;
for j=2:(length(pts)-1)
    tmpcurve = ReSampleCurve(pts(:,j:j+1),closest(j+1)-closest(j)+1);
    intcurve = [intcurve,tmpcurve(:,2:end)];
end
tmpcurve = ReSampleCurve(pts(:,[k,1]),N-closest(k)+closest(1));
intcurve = [intcurve,tmpcurve(:,2:end)];

% Shifts coordinates of original curve to match first point of linear
% reconstruction
if (closest(1)==1)
    Xn = X(:,:);
else
    Xn = ShiftF(X(:,1:(N-1)),closest(1)-1);
    Xn(:,end+1) = Xn(:,1);
end

% Compute linear reconstruction error
qi = curve_to_q(intcurve);
q = curve_to_q(Xn);
dst = (norm(q(1,:)-qi(1,:))^2)+(norm(q(2,:)-qi(2,:))^2);

% Display original curve and linear reconstruction
if fig==1
    hold on
    axis off
    plot(X(1,:),X(2,:),'LineWidth',3)
    plot(intcurve(1,:),intcurve(2,:),'g','LineWidth',3)
    scatter(pts(1,:),pts(2,:),100,'Filled')
end