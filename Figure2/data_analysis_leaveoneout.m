% Perform no-leave-out analysis on the cluster
%clear all

AI=getenv('PBS_ARRAYID');
job_dex=str2num(AI);
intervals=[[1 3];[4 6];[7 9];[10 12]];
st=intervals(job_dex,1);
stend=intervals(job_dex,2);

r2dataeb=0; r2datashuf=1; r2datasurr=0;
fprintf('Analyzing upper bound data...')
solve_leaveoneout_fxn(st,stend,r2dataeb,r2datashuf,r2datasurr);
%save -v7.3 leaveoneout_surr
%prune_loo(0,0,1);

exit
