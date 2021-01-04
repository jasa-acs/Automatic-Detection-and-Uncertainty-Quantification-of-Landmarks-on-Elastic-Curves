function q = curve_to_q(p)
%% Converts a curve p to its SRVF representation q
%% Input
% p = one two-dimensional curve
%% Outputs
% q = SRVF of p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute time-derivative of p
[d,N] = size(p);
for i = 1:d
    v(i,:) = gradient(p(i,:),1/(N-1));
end

% Computes SRVF of p
for i = 1:N
    L(i) = sqrt(norm(v(:,i)));
    if L(i) > 0.0001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = v(:,i)*0.0001;
    end
end