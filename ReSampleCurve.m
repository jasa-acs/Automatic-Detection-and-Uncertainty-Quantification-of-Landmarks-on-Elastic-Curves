function Xn = ReSampleCurve(X,N)
%% Resamples X to N points (arc-length parameterized)
%% Inputs
% X = two-dimensional curve
% N = # of discretization points desired
%% Outputs
% Xn = resampled curve of X with N points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,T] = size(X);

% Computes distance between neighboring points (spacing) on X
del(1) = 0;
for r = 2:T
    del(r) = norm(X(:,r) - X(:,r-1));
end
cumdel = cumsum(del)/sum(del);

% Desired new spacing (with N points)
newdel = [0:(N-1)]/(N-1);

% Interpolate X at new values newdel
for j=1:d
    Xn(j,:) = interp1(cumdel,X(j,1:T),newdel,'linear');
end
