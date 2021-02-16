% Perform no-leave-out data analysis from recording session st to recording
% session stend
function solve_noleaveout_fxn(st,stend)

load covvarcorrdata

p = pwd;
addpath(p)
addpath(fullfile(p,'lib'))
addpath(fullfile(p,'interface'))

r2dataeb=1; % value of 1 means original data will be analyzed
r2datashuf=0; % value of 1 means shuffled data will be analyzed
r2datasurr=0; % value of 1 means upper bound data will be analyzed
nshuf=100; % number of shuffles to be performed
mu=81; % maximum number of neurons
nr=100; % number of Poisson processes to be generated for upper bound
ndata=37; % number of recording sessions

% save all covariance matrices you will use: attended, unattended, actual
% data, shuffled data, upper bound data, right hemisphere, left hemisphere
covmatallrealUR=zeros(ndata,mu,mu)+NaN;
covmatallshufUR=zeros(ndata,nshuf,mu,mu)+NaN;
covmatallsurrUR=zeros(ndata,nr,mu,mu)+NaN;
covmatallrealAR=zeros(ndata,mu,mu)+NaN;
covmatallshufAR=zeros(ndata,nshuf,mu,mu)+NaN;
covmatallsurrAR=zeros(ndata,nr,mu,mu)+NaN;
covmatallrealUL=zeros(ndata,mu,mu)+NaN;
covmatallshufUL=zeros(ndata,nshuf,mu,mu)+NaN;
covmatallsurrUL=zeros(ndata,nr,mu,mu)+NaN;
covmatallrealAL=zeros(ndata,mu,mu)+NaN;
covmatallshufAL=zeros(ndata,nshuf,mu,mu)+NaN;
covmatallsurrAL=zeros(ndata,nr,mu,mu)+NaN;

nunknownmax=mu; % maximum number of unknowns to be solved for
nknownmax=mu*(mu-1)/2; % maximum number of covariance values

% save all g-values for actual data "real", shuffled data "shuf", upper
% bound data "surr", for right and left hemispheres, for each recording
% session
xvecallrealdataR=zeros(ndata,nunknownmax)+NaN;
xvecallshufdataR=zeros(ndata,nshuf,nunknownmax)+NaN;
xvecallsurrdataR=zeros(ndata,nr,nunknownmax)+NaN;
xvecallrealdataL=zeros(ndata,nunknownmax)+NaN;
xvecallshufdataL=zeros(ndata,nshuf,nunknownmax)+NaN;
xvecallsurrdataL=zeros(ndata,nr,nunknownmax)+NaN;

% These matrices contain the information for the sets of equations that are
% being solved. First dimension refers to the recording session, second
% dimension refers to either which equation (one for each pair of neurons),
% or which shuffle or Poisson process we are on, and the last dimension
% contains the four values that define each equation (explained later)
eqindvecrealR=zeros(ndata,nknownmax,4)+NaN;
eqindvecshufR=zeros(ndata,nshuf,nknownmax,4)+NaN;
eqindvecsurrR=zeros(ndata,nr,nknownmax,4)+NaN;
eqindvecrealL=zeros(ndata,nknownmax,4)+NaN;
eqindvecshufL=zeros(ndata,nshuf,nknownmax,4)+NaN;
eqindvecsurrL=zeros(ndata,nr,nknownmax,4)+NaN;

% Save numbers of covariance values used for each recording session
ncomprealR=zeros(37,1);
ncompshufR=zeros(37,1);
ncompsurrR=zeros(37,1);
ncomprealL=zeros(37,1);
ncompshufL=zeros(37,1);
ncompsurrL=zeros(37,1);

% Save r^2 values for each recording session
r2realR=zeros(37,1);
r2shufR=zeros(37,1);
r2surrR=zeros(37,1);
r2realL=zeros(37,1);
r2shufL=zeros(37,1);
r2surrL=zeros(37,1);

% Iterate through all recording sessions (or just one of them, if you
% un-comment "for ndi=10" and comment "for ndi=st:stend")
%st=1;
%stend=ndata;
for ndi=st:stend
%for ndi=10
    comprealR=[]; % These will contain cov values (actual vs predicted)
    compshufR=[];
    compsurrR=[];
    comprealL=[];
    compshufL=[];
    compsurrL=[];
    fprintf('Analyzing data set %d out of %d ... \n',ndi,ndata);
    % Right hemisphere
    nru=sum(~isnan(allcovmatunR(:,1,ndi))); % number of units in right hemisphere for this rec. session
    % Extract cov matrices for this rec. session
    covmatunatt=squeeze(allcovmatunR(1:nru,1:nru,ndi));
    covmatatt=squeeze(allcovmatatR(1:nru,1:nru,ndi));
    % save exact cov matrix used for this analysis
    covmatallrealUR(1:nru,1:nru,ndi)=covmatunatt;
    covmatallrealAR(1:nru,1:nru,ndi)=covmatatt;
    % temporarily save original cov matrices used (before they are shuffled
    % etc)
    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    if r2dataeb==1 % analyze actual data
        fprintf('Computing r2 of real data ... \n')
        tic
        solvecovmat; % Equations solved in here
        xvecallrealdataR(ndi,1:nru)=g; % Save all g-values
        eqindvecrealR(ndi,1:ncv,:)=eqind; % Save system of equations
        compare_r2; % Get actual vs predicted cov values
        toc
        comprealR=[comprealR;[cprevec cvalvec]]; % Save actual vs predicted cov values
    end
    if r2datashuf==1 % analyze shuffled data
        fprintf('Computing r2 of shuffled data ... \n')
        for nsi=1:nshuf % do nshuf different shuffles
            % for each of them, start out with original cov matrices
            covmatunatt=covmatunatt_parent;  
            covmatatt=covmatatt_parent;
            %tic
            shuffle_covmat; % shuffle the cov matrices
            solvecovmat; % solve for g
            covmatallshufUR(ndi,nsi,1:nru,1:nru)=covmatunatt;
            covmatallshufAR(ndi,nsi,1:nru,1:nru)=covmatatt;
            eqindvecshufR(ndi,nsi,1:ncv,:)=eqind;
            compare_r2;
            %toc
            xvecallshufdataR(ndi,nsi,1:nru)=g;
            compshufR=[compshufR;[cprevec cvalvec]];
        end
    end
    if r2datasurr==1 % analyze upper bound
        fprintf('Computing r2 of surrogate data set %d out of %d ... \n',ndi,ndata)
        covmatunatt=covmatunatt_parent; % start with original cov matrices
        covmatatt=covmatatt_parent;
        solvecovmat; % solve for g
        covmatatt_p=zeros(nru);
        % construct covmatatt_p, for which g provides an exact solution
        for f1=1:nru
            for f2=1:nru
                covmatatt_p(f1,f2)=g(f1)*g(f2)*covmatunatt(f1,f2);
            end
        end
        % variance values
        muat=diag(covmatatt_p);
        % number of samples
        Nsamples=ngoodtratR(ndi);    
        % Prepare to generate correlated Poisson processes (only need to
        % run once to save time)
        [pmfs,supports] = PoissonMarginals(muat,1e-8);
        [gammas,Lambda,joints2D]=FindDGAnyMarginal(pmfs,covmatatt_p,supports);
        for ri=1:nr % generate nr separate Poisson processes
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            if rem(ri,20)==0
                fprintf('Run number %d for data set %d \n',ri,ndi)
            end
            % Generate Poisson processes
            [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples);
            if length(Sat)==1
                break; % Uh oh, didn't work. Probably a zero somewhere
            end
            Sat = Sat';
            muathat=mean(Sat,2);
            Cathat=cov(Sat');
            covmatatt=Cathat;
            solvecovmat; % Solve for g
            compare_r2;
            covmatallsurrUR(ndi,ri,1:nru,1:nru)=covmatunatt;
            covmatallsurrAR(ndi,ri,1:nru,1:nru)=covmatatt;
            eqindvecsurrR(ndi,ri,1:ncv,:)=eqind;
            xvecallsurrdataR(ndi,ri,1:nru)=g;
            compsurrR=[compsurrR;[cprevec cvalvec]];
        end
    end
    
    %Left Hemisphere
    nlu=sum(~isnan(allcovmatunL(:,1,ndi)));
    covmatunatt=squeeze(allcovmatunL(1:nlu,1:nlu,ndi));
    covmatatt=squeeze(allcovmatatL(1:nlu,1:nlu,ndi));
    covmatallrealUL(1:nlu,1:nlu,ndi)=covmatunatt;
    covmatallrealAL(1:nlu,1:nlu,ndi)=covmatatt;
    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    if r2dataeb==1
        fprintf('Left hemisphere ... \n')
        solvecovmat;
        xvecallrealdataL(ndi,1:nlu)=g;
        eqindvecrealL(ndi,1:ncv,:)=eqind;
        compare_r2;
        comprealL=[comprealL;[cprevec cvalvec]];
    end
    if r2datashuf==1
        fprintf('Left hemisphere... \n')
        for nsi=1:nshuf
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            shuffle_covmat;
            solvecovmat;
            covmatallshufUL(ndi,nsi,1:nlu,1:nlu)=covmatunatt;
            covmatallshufAL(ndi,nsi,1:nlu,1:nlu)=covmatatt;
            eqindvecshufL(ndi,nsi,1:ncv,:)=eqind;
            compare_r2;
            xvecallshufdataL(ndi,nsi,1:nlu)=g;
            compshufL=[compshufL;[cprevec cvalvec]];
        end
    end
    if r2datasurr==1
        fprintf('Left hemisphere... \n')
        covmatunatt=covmatunatt_parent;
        covmatatt=covmatatt_parent;
        solvecovmat;
        covmatatt_p=zeros(nlu);
        for f1=1:nlu
            for f2=1:nlu
                covmatatt_p(f1,f2)=g(f1)*g(f2)*covmatunatt(f1,f2);
            end
        end
        muat=diag(covmatatt_p);
        Nsamples=ngoodtratL(ndi);           
        [pmfs,supports] = PoissonMarginals(muat,1e-8);
        [gammas,Lambda,joints2D]=FindDGAnyMarginal(pmfs,covmatatt_p,supports);
        for ri=1:nr
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            if rem(ri,20)==0
                fprintf('Run number %d for data set %d \n',ri,ndi)
            end
            %tic
            [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples);
            if length(Sat)==1
                break;
            end
            Sat = Sat';
            muathat=mean(Sat,2);
            Cathat=cov(Sat');
            covmatatt=Cathat;
            solvecovmat;
            compare_r2;
            covmatallsurrUL(ndi,ri,1:nlu,1:nlu)=covmatunatt;
            covmatallsurrAL(ndi,ri,1:nlu,1:nlu)=covmatatt;
            eqindvecsurrL(ndi,ri,1:ncv,:)=eqind;
            xvecallsurrdataL(ndi,ri,1:nlu)=g;
            compsurrL=[compsurrL;[cprevec cvalvec]];
        end
    end 
    if r2dataeb
        r2rR=corrcoef(comprealR);
        r2realR(ndi)=r2rR(1,2);
        fprintf('r2 of RH real data set %d is %f \n',ndi,r2realR(ndi))
        ncomprealR(ndi)=length(comprealR);
        r2rL=corrcoef(comprealL);
        r2realL(ndi)=r2rL(1,2);
        fprintf('r2 of LH real data set %d is %f \n',ndi,r2realL(ndi))
        ncomprealL(ndi)=length(comprealL);
    end
    if r2datashuf
        r2shR=corrcoef(compshufR);
        r2shufR(ndi)=r2shR(1,2);
        fprintf('r2 of RH shuffled data set %d is: %f \n',ndi,r2shufR(ndi))
        ncompshufR(ndi)=length(compshufR);
        r2shL=corrcoef(compshufL);
        r2shufL(ndi)=r2shL(1,2);
        fprintf('r2 of LH shuffled data set %d is: %f \n',ndi,r2shufL(ndi))
        ncompshufL(ndi)=length(compshufL);
    end
    if r2datasurr
        r2srR=corrcoef(compsurrR);
        r2surrR(ndi)=r2srR(1,2);
        fprintf('r2 of RH surr data set %d is %f \n',ndi,r2surrR(ndi))
        ncompsurrR(ndi)=length(compsurrR);
        r2srL=corrcoef(compsurrL);
        r2surrL(ndi)=r2srL(1,2);
        fprintf('r2 of LH surr data set %d is %f \n',ndi,r2surrL(ndi))
        ncompsurrL(ndi)=length(compsurrL);
    end

    %sname=['nlo_nunits' num2str(ndi) '.mat'];
    %save(sname);
    save -v7.3 sname
    %if ndi>st
    %  dname=['nlo_nunits' num2str(ndi-1) '.mat'];
    %  delete(dname)
    %end
end

end
