clear all

p = pwd;
addpath(p)
addpath(fullfile(p,'lib'))
addpath(fullfile(p,'interface'))

load covvarcorrdata

sampvec=10:100:10000;
perfvec=zeros(2,ndata);
extrvec=zeros(2,ndata,length(sampvec));
for ndi=1:ndata
    ndi
    nru=sum(squeeze(~isnan(allcovmatunR(1,:,ndi))));
    covmatunatt=squeeze(allcovmatunR(1:nru,1:nru,ndi));
    covmatatt=squeeze(allcovmatatR(1:nru,1:nru,ndi));

    covmatunatt_parent=covmatunatt;
    covmatatt_parent=covmatatt;
    
    solvecovmat;
    covmatatt_p=zeros(nru);
    for f1=1:nru
        for f2=1:nru
            covmatatt_p(f1,f2)=g(f1)*g(f2)*covmatunatt(f1,f2);
        end
    end
    muat=diag(covmatatt_p);
    Nsamples=ngoodtratR(ndi);
    [pmfs,supports] = PoissonMarginals(muat,1e-8);
    [gammas,Lambda,joints2D]=FindDGAnyMarginal(pmfs,covmatatt_p,supports);
    [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples);
    Cathat=cov(Sat);
    covmatatt=Cathat;
    solvecovmat;
    compare_r2;
    perfvec(1,ndi)=Nsamples;
    perfvec(2,ndi)=r2best(1,2);
    nnni=1;
    for nnn=sampvec
        nnn
        [Sat,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,nnn);
        Cathat=cov(Sat);
        covmatatt=Cathat;
        solvecovmat;
        compare_r2;
        extrvec(1,ndi,nnni)=nnn;
        extrvec(2,ndi,nnni)=r2best(1,2);
        nnni=nnni+1;
    end
end

save surrconvdata