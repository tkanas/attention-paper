%% Leave-one-out analysis
function solve_leaveoneout1_fxn(st,stend,r2dataeb,r2datashuf,r2datasurr)

load covvarcorrdata

p = pwd;
addpath(p)
addpath(fullfile(p,'lib'))
addpath(fullfile(p,'interface'))

% Run each of these options separately so as not to run out of memory
%r2dataeb=0; % actual data
%r2datashuf=0; % shuffled data
%r2datasurr=1; % upper bound data

nshuf=100; % number of shuffles
mu=81; % max number of neurons
nr=100; % number of Poisson processes
ndraws=1000; % number of leave-one-outs

r2realR=zeros(ndata,1);
r2shufR=zeros(ndata,nshuf);
r2surrR=zeros(ndata,nr);
r2realL=zeros(ndata,1);
r2shufL=zeros(ndata,nshuf);
r2surrL=zeros(ndata,nr);

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

fprintf('Entered solve routine \n')

nunknowmu=mu;
if r2dataeb>0
    xvecallrealdataR=zeros(ndata,ndraws,nunknowmu)+NaN;
    eqremvecrealR=zeros(ndata,ndraws,4)+NaN;
    xvecallrealdataL=zeros(ndata,ndraws,nunknowmu)+NaN;
    eqremvecrealL=zeros(ndata,ndraws,4)+NaN;
    comprealRall=cell(37,1);
    comprealLall=cell(37,1);
end
if r2datashuf>0
    xvecallshufdataR=zeros(ndata,nshuf,ndraws,nunknowmu)+NaN;
    eqremvecshufR=zeros(ndata,nshuf,ndraws,4)+NaN;
    xvecallshufdataL=zeros(ndata,nshuf,ndraws,nunknowmu)+NaN;
    eqremvecshufL=zeros(ndata,nshuf,ndraws,4)+NaN;
    compshufRall=cell(37,1);
    compshufLall=cell(37,1);
end
if r2datasurr>0
    xvecallsurrdataR=zeros(ndata,nr,ndraws,nunknowmu)+NaN;
    eqremvecsurrR=zeros(ndata,nr,ndraws,4)+NaN;
    xvecallsurrdataL=zeros(ndata,nr,ndraws,nunknowmu)+NaN;
    eqremvecsurrL=zeros(ndata,nr,ndraws,4)+NaN;
    compsurrRall=cell(37,1);
    compsurrLall=cell(37,1);
end

ncomprealR=zeros(37,1);
ncompshufR=zeros(37,1);
ncompsurrR=zeros(37,nr);
ncomprealL=zeros(37,1);
ncompshufL=zeros(37,1);
ncompsurrL=zeros(37,nr);


for ndi=st:stend
%for ndi=1:6
    comprealR=[];
    compshufR=[];
    compsurrR=[];
    comprealL=[];
    compshufL=[];
    compsurrL=[];
    fprintf('Analyzing data set %d out of %d ... \n',ndi,ndata);
    
    % Right hemisphere
    nru=sum(~isnan(allcovmatunR(:,1,ndi)));

    covmatunatt=squeeze(allcovmatunR(1:nru,1:nru,ndi));
    covmatatt=squeeze(allcovmatatR(1:nru,1:nru,ndi));
    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    ndraws=min(1000,nru*(nru-1)/2);
    if r2dataeb==1
        fprintf('Computing r2 of RH real data ... \n')
        covmatallrealUR(1:nru,1:nru,ndi)=covmatunatt;
        covmatallrealAR(1:nru,1:nru,ndi)=covmatatt;
        solvecovmat_loo; % Solve system of equations
        eqremvecrealR(ndi,1:ndraws,:)=eqremvec;
        compare_r2_loo;
        xvecallrealdataR(ndi,1:ndraws,1:nru)=gvecall;
        comprealR=[comprealR;[cprevec cvalvec]];
        comprealRall{ndi}=comprealR;
        r2rR=corrcoef(comprealR);
        r2realR(ndi)=r2rR(1,2);
        ncomprealR(ndi)=length(comprealR);
    end
    if r2datashuf==1
        fprintf('Computing r2 of RH shuffled data ... \n')
        cmpshfRall=[];
        for nsi=1:nshuf
            compshufR=[];
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            shuffle_covmat;
            solvecovmat_loo; % solve system of equations
            covmatallshufUR(ndi,nsi,1:nru,1:nru)=covmatunatt;
            covmatallshufAR(ndi,nsi,1:nru,1:nru)=covmatatt;
            eqremvecshufR(ndi,nsi,1:ndraws,:)=eqremvec;
            compare_r2_loo; 
            xvecallshufdataR(ndi,nsi,1:ndraws,1:nru)=gvecall;
            compshufR=[compshufR;[cprevec cvalvec]];
            r2shR=corrcoef(compshufR);
            r2shufR(ndi,nsi)=r2shR(1,2);
            cmpshfRall=[cmpshfRall;[cprevec cvalvec]];
        end
        compshufRall{ndi}=cmpshfRall;
    end
    if r2datasurr==1
        fprintf('Computing r2 of RH surrogate data set %d out of %d ... \n',ndi,ndata)
        covmatunatt=covmatunatt_parent;
        covmatatt=covmatatt_parent;
        solvecovmat;
        covmatatt_p=zeros(nru);
        for f1=1:nru
            for f2=1:nru
                covmatatt_p(f1,f2)=g(f1)*g(f2)*covmatunatt(f1,f2);
            end
        end
        muat=diag(covmatatt_p);
        Nsamples=ngoodtratR(ndi);
        fprintf('Generating info of DG... \n')
        [pmfs,supports] = PoissonMarginals(muat,1e-8);
        [gammas,Lambda,joints2D]=FindDGAnyMarginal(pmfs,covmatatt_p,supports);
        cmpsurrRall=[];
        fprintf('Starting different Poiss P... \n')
        for ri=1:nr
            compsurrR=[];
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            if rem(ri,20)==0
                fprintf('Run number %d for data set %d \n',ri,ndi)
            end
            fprintf('Maybe the problem ... \n')
            [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples);
            if length(Sat)==1
                break;
            end
            Sat=Sat';
            muathat=mean(Sat,2);
            Cathat=cov(Sat');
            covmatatt=Cathat;
            solvecovmat_loo;
            compare_r2_loo;
            covmatallsurrUR(ndi,ri,1:nru,1:nru)=covmatunatt;
            covmatallsurrAR(ndi,ri,1:nru,1:nru)=covmatatt;
            eqremvecsurrR(ndi,ri,1:ndraws,:)=eqremvec;
            xvecallsurrdataR(ndi,ri,1:ndraws,1:nru)=gvecall;
            cmpsurrRall=[cmpsurrRall;[cprevec cvalvec]];
            compsurrR=[compsurrR;[cprevec cvalvec]];
            r2srR=corrcoef(compsurrR);
            r2surrR(ndi,ri)=r2srR(1,2);
            ncompsurrR(ndi,ri)=length(compsurrR);
        end
        compsurrRall{ndi}=cmpsurrRall;
    end
     
    % Left hemisphere
    nlu=sum(~isnan(allcovmatunL(:,1,ndi)));
    
    covmatunatt=squeeze(allcovmatunL(1:nlu,1:nlu,ndi));
    covmatatt=squeeze(allcovmatatL(1:nlu,1:nlu,ndi));
    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    ndraws=min(1000,nlu*(nlu-1)/2);
    if r2dataeb==1
        comprealL=[];
        fprintf('Computing r2 of LH real data ... \n')
        covmatallrealUL(1:nlu,1:nlu,ndi)=covmatunatt;
        covmatallrealAL(1:nlu,1:nlu,ndi)=covmatatt;
        solvecovmat_loo;
        eqremvecrealL(ndi,1:ndraws,:)=eqremvec;
        compare_r2_loo;
        xvecallrealdataL(ndi,1:ndraws,1:nlu)=gvecall;
        comprealL=[comprealL;[cprevec cvalvec]];
        comprealLall{ndi}=comprealL;
        r2rL=corrcoef(comprealL);
        r2realL(ndi)=r2rL(1,2);
        ncomprealL(ndi)=length(comprealL);
    end
    if r2datashuf==1
        fprintf('Computing r2 of LH shuffled data ... \n')
        cmpshfLall=[];
        for nsi=1:nshuf
            compshufL=[];
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
            shuffle_covmat;
            solvecovmat_loo;
            covmatallshufUL(ndi,nsi,1:nlu,1:nlu)=covmatunatt;
            covmatallshufAL(ndi,nsi,1:nlu,1:nlu)=covmatatt;
            eqremvecshufL(ndi,nsi,1:ndraws,:)=eqremvec;
            compare_r2_loo;
            xvecallshufdataL(ndi,nsi,1:ndraws,1:nlu)=gvecall;
            compshufL=[compshufL;[cprevec cvalvec]];
            r2shL=corrcoef(compshufL);
            r2shufL(ndi,nsi)=r2shL(1,2);
            cmpshfLall=[cmpshfLall;[cprevec cvalvec]];
        end
        compshufLall{ndi}=cmpshfLall;
    end
    if r2datasurr==1
        fprintf('Computing r2 of LH surrogate data set %d out of %d ... \n',ndi,ndata)
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
        cmpsurrLall=[];
        for ri=1:nr
            compsurrL=[];
            covmatunatt=covmatunatt_parent;
            covmatatt=covmatatt_parent;
           % if rem(ri,20)==0
           %     fprintf('Run number %d for data set %d \n',ri,ndi)
           % end
            [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples);
            if length(Sat)==1
                break;
            end
            Sat=Sat';
            muathat=mean(Sat,2);
            Cathat=cov(Sat');
            covmatatt=Cathat;
            solvecovmat_loo;
            compare_r2_loo;
            covmatallsurrUL(ndi,ri,1:nlu,1:nlu)=covmatunatt;
            covmatallsurrAL(ndi,ri,1:nlu,1:nlu)=covmatatt;
            eqremvecsurrL(ndi,ri,1:ndraws,:)=eqremvec;
            xvecallsurrdataL(ndi,ri,1:ndraws,1:nlu)=gvecall;
            compsurrL=[compsurrL;[cprevec cvalvec]];
            cmpsurrLall=[cmpsurrLall;[cprevec cvalvec]];
            r2srL=corrcoef(compsurrL);
            r2surrL(ndi,ri)=r2srL(1,2);
            ncompsurrL(ndi,ri)=length(compsurrL);
        end
        compsurrLall{ndi}=cmpsurrLall;
    end
    %{
    if r2dataeb==1
        r2rR=corrcoef(comprealR);
        r2realR(ndi)=r2rR(1,2);
        fprintf('r2 of RH real data set %d is %f \n',ndi,r2realR(ndi))
        ncomprealR(ndi)=length(comprealR);
        r2rL=corrcoef(comprealL);
        r2realL(ndi)=r2rL(1,2);
        fprintf('r2 of LH real data set %d is %f \n',ndi,r2realL(ndi))
        ncomprealL(ndi)=length(comprealL);
    end
    if r2datashuf==1
        r2shR=corrcoef(compshufR);
        r2shufR(ndi)=r2shR(1,2);
        fprintf('r2 of RH shuffled data set %d is: %f \n',ndi,r2shufR(ndi))
        ncompshufR(ndi)=length(compshufR);
        r2shL=corrcoef(compshufL);
        r2shufL(ndi)=r2shL(1,2);
        fprintf('r2 of LH shuffled data set %d is: %f \n',ndi,r2shufL(ndi))
        ncompshufL(ndi)=length(compshufL);
    end
    
    if r2datasurr==1
        r2srR=corrcoef(compsurrR);
        r2surrR(ndi)=r2srR(1,2);
        fprintf('r2 of RH surr data set %d is %f \n',ndi,r2surrR(ndi))
        ncompsurrR(ndi)=length(compsurrR);
        r2srL=corrcoef(compsurrL);
        r2surrL(ndi)=r2srL(1,2);
        fprintf('r2 of LH surr data set %d is %f \n',ndi,r2surrL(ndi))
        ncompsurrL(ndi)=length(compsurrL);
    end
    %}
    sname=['loo_units' num2str(ndi) '.mat'];
    save(sname,'-v7.3');
    if ndi>st
      dname=['loo_units' num2str(ndi-1) '.mat'];
      delete(dname);
    end
end

end
