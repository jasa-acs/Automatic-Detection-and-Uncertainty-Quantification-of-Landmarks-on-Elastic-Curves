function pn = ShiftF(p,tau)
%% Shifts coordinate ordering of p
%% Inputs
% p = two-dimensional curve
% tau = shift amount (to the left)
%% Outputs
% pn = new shifted coordinates of p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,N] = size(p);

% No shift
if (tau == 0)
    pn = p;
    return;
end

% Nonzero shift
if (tau > 0)
    pn(:,1:(N-tau)) = p(:,(tau+1):N);
    pn(:,(N-tau+1):N) = p(:,1:tau);
else
    t = abs(tau)+1;
    pn(:,1:(N-t+1)) = p(:,t:N);
    pn(:,(N-t+2):N) = p(:,1:(t-1));
end