% Prune leave-one-out data.
% Input: triplet with one 1 and two 0's, indicating which data to prune.
% Must update the "load" commands to whatever the data is called
function [] = prune_loo(prune_real,prune_shuf,prune_surr)
    if prune_real==1
        filelist=dir('loo_nunits_real*.mat');
        YAmR=zeros(37,1);
        YAmL=zeros(37,1);
        a1=zeros(37,1);
        for f=1:length(filelist)
            load(filelist(f).name);
            for i=st:stend
                YAmR(i)=r2realR(i);
                YAmL(i)=r2realL(i);
                a1(i)=i;
            end
        end
        save pruned_loo_real a1 YAmR YAmL
    elseif prune_shuf==1
        filelist=dir('loo_nunits_shuf*.mat');
        YSmR=zeros(37,1); YSseR=zeros(37,1);
        YSmL=zeros(37,1); YSseL=zeros(37,1);
        a1=zeros(37,1);
        for f=1:length(filelist)
            load(filelist(f).name);
            for i=st:stend
                vR=r2shufR(i,:);
                YSmR(i)=mean(vR);
                YSseR(i)=std(vR)/sqrt(length(vR));
                vL=r2shufL(i,:);
                YSmL(i)=mean(vL);
                YSseL(i)=std(vL)/sqrt(length(vL));
                a1(i)=i;
            end
        end
        save pruned_loo_shuf a1 YSmR YSseR YSmL YSseL
    elseif prune_surr==1
        filelist=dir('loo_nunits_surr*.mat');
        YUmR=zeros(37,1); YUseR=zeros(37,1);
        YUmL=zeros(37,1); YUseL=zeros(37,1);
        a1=zeros(37,1);
        for f=1:length(filelist)
            load(filelist(f).name);
            for i=st:stend
                vR=r2surrR(i,:);
                YUmR(i)=mean(vR);
                YUseR(i)=std(vR)/sqrt(length(vR));
                vL=r2surrL(i,:);
                YUmL(i)=mean(vL);
                YUseL(i)=std(vL)/sqrt(length(vL));
                a1(i)=i;
            end
        end
        save pruned_loo_surr a1 YUmR YUseR YUmL YUseL
    end
end