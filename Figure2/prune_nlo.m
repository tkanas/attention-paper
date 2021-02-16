% Prune dataset to save only values to be plotted. This is important for
% larger datasets like for leave-one-out and partition. 

clear all

filelist=dir('nlo_nunits*.mat');

YSmR=[]; YSseR=[];
YSmL=[]; YSseL=[];
YUmR=[]; YUseR=[];
YUmL=[]; YUseL=[];

for f=1:length(filelist)
    load(filelist(f).name);
    for i=st:stend
        YSmR=[YSmR;mean(r2shufR(i,:))];
        YSseR=[YSseR;std(r2shufR(i,:))/sqrt(length(r2shufR(i,:)))];
        YSmL=[YSmL;mean(r2shufL(i,:))];
        YSseL=[YSseL;std(r2shufL(i,:))/sqrt(length(r2shufL(i,:)))];
        YUmR=[YUmR;mean(r2surrR(i,:))];
        YUseR=[YUseR;std(r2surrR(i,:))/sqrt(length(r2surrR(i,:)))];
        YUmL=[YUmL;mean(r2surrL(i,:))];
        YUseL=[YUseL;std(r2surrL(i,:))/sqrt(length(r2surrL(i,:)))];
    end
end
        
save noleaveout_all YSmR YSseR YSmL YSseL YUmR YUseR YUmL YUseL r2realR r2realL