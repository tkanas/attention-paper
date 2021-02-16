% Perform no-leave-out analysis on the cluster
clear all

%load covvarcorrdata
%load covmats_lr_new

% No-leave-out analysis
%tic

AI=getenv('PBS_ARRAYID');
job_dex=str2num(AI);
intervals=[[1 4];[5 7];[8 10];[11 13];[14 16];[17 19];[20 22];[23 25];[26 28];[29 31];[32 34];[35 37]];
st=intervals(job_dex,1);
stend=intervals(job_dex,2);

solve_noleaveout_fxn(st,stend);
%save noleaveout
%toc
%save -v7.3 noleaveout
exit
