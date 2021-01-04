function val = InnerProd_Q(q1,q2)
%% L2 inner product of q1 and q2
%% Input
% q1, q2 = two-dimensional curves
%% Output
% val = inner product of q1 and q2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,N] = size(q1);
val = trapz(linspace(0,1,N),sum(q1.*q2));