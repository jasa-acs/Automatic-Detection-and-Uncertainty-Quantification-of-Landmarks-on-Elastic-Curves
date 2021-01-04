%% How to use the included code
%% Section 5 (simulated data)
load('ToyCurve.mat')    % loads simulated data - description found in ToyCurvelabels.txt

% In order to reproduce 5.1, use the following (with default settings) -
% description of how to alter hyperparameters/settings found in ALDfixed.m
% comments
ALDfixed(X200,4,1,1)

% 5.2 can be performed by simply rotating and re-scaling X200 according to
% the manuscript, and then repeating the above
theta = pi/4;
R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
RotX200 = R*X200;
ScaX200 = 2*X200;

ALDfixed(RotX200,4,1,1)
ALDfixed(ScaX200,4,1,1)

% The simulation study of 5.3 simply requires changing the hyperparameters
% a and b to reflect those of the table, e.g.,
ALDfixed(X200,4,1,1,5,0.01,[],[],[])

% Section 5.4 uses the FindBestk function with parallel computing and
% compares it to RJMCMC
FindBestk(X100,10,[],[],[],1e5,1e6,[])  % uses d_k^2 criterion - description of settings found in comments of FindBestk.m

ALDunknown(X100,1,1,1)      % uses RJMCMC with lambda = 1
ALDunknown(X100,1e-6,1,1)   % lambda = 10^-6

% Section 5.5 uses the set of 5 similar simulated curves to perform fixed k
% inference
ALDfixed(MultX200,4,1,1)

%% Section 6 (real-life examples) -- this is more free-flowing, but some example lines to execute are included below
% Bone from MPEG-7 in Section 6.1
load('MPEG7closed.mat')
ALDfixed(C_cl(:,:,1),4,1,1)     % only alter the last argument of C_cl to change the particular shape

% All 20 bones from MPEG-7
ALDfixed(C_cl(:,:,1:20),4,1,1)

% RJMCMC with deer from MPEG-7
ALDunknown(C_cl(:,:,461),0.00001,1,1)       % lambda = 0.00001
ALDunknown(C_cl(:,:,461),0.001,1,1)         % lambda = 0.001

% All 30 mice vertebrae from Control group
load('mice.mat')
ALDfixed(MO(:,:,1:30),4,1,1)

%% Note: If any of the inputs are unclear, the corresponding .m file should have detailed comments describing the inputs and outputs generated