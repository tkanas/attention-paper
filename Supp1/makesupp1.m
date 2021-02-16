clear all
load covvarcorrdata

Nvec=10:10:150;
Nvec=[Nvec 500];
numsamp=length(Nvec);
r2vec=zeros(ndata,numsamp);
nunits=zeros(ndata,1);

ns=100; % number of shuffles

for ndi=1:ndata
    ndi
    nru=sum(squeeze(~isnan(allcovmatatR(1,:,ndi))));
    nunits(ndi)=nru;
    covmatunatt=squeeze(allcovmatunR(1:nru,1:nru,ndi));
    covmatatt=squeeze(allcovmatatR(1:nru,1:nru,ndi));
    
    % Save original cov matrices
    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    
    % Turn cov matrices into vectors
    cmuv=reshape(covmatunatt,numel(covmatunatt),1);
    cmav=reshape(covmatatt,numel(covmatatt),1);
    
    % Mean and std of cov vectors
    mcunat=mean(cmuv);
    mcat=mean(cmav);
    stdunat=std(cmuv);
    stdat=std(cmav);
    
    % Generate synthetic data (random Gaussian) for various numbers of
    % units
    nni=1;
    for N=Nvec
        N
        comparevec=[];
        for nsi=1:ns
            covmatunatt=randn(N)*stdunat+mcunat;
            covmatatt=randn(N)*stdat+mcat;
            solvecovmat;
            compare_r2;
            comparevec=[comparevec;[cvalvec cprevec]];
        end
        r2b=corrcoef(comparevec);
        r2vec(ndi,nni)=r2b(1,2);
        nni=nni+1;
    end
end
save shufconvdata1