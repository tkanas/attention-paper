% Run the no-leave-out analysis NOT on the cluster (e.g. if computing
% single example)
clear all

load covvarcorrdata

% No-leave-out analysis
tic
solve_noleaveout;
save noleaveout
toc