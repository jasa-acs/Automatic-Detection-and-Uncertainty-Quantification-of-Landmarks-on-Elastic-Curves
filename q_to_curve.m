function p = q_to_curve(q)
%% Converts a SRVF q to its coordinate function p
%% Input
% q = one two-dimensional SRVF
%% Outputs
% p = coordinate function of q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes Euclidean norm of SRVF q as function of curve domain
[d,N] = size(q);
for i = 1:N
    qnorm(i) = norm(q(:,i),'fro');
end

% Computes coordinate function of q
for i = 1:d
    p(i,:) = [cumtrapz(linspace(0,1,N),q(i,:).*qnorm)];
end