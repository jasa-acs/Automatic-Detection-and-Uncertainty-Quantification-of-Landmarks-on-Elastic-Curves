function p = MoveProb(k,k_min)
%% Computes probability of birth/death/stay at current value of k
%% Inputs
% Required:
% k = current number of landmarks
% k_min = minimum number of landmarks allowed
%% Output
% p = birth/death/stay probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Birth
p_b = 1/3;

% Death
if k > k_min
    p_d = 1/3;
else
    p_d = 0;
end

% Stay
p_s = 1-p_b-p_d;

p = [p_b;p_d;p_s];