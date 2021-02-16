% Plot Figure 2 h, i
% Bar plots for leave-one-out cross-validation
clear all

% Pruned leave-one-out data
load pruned_loo_real
load pruned_loo_shuf
load pruned_loo_surr

% Bar color scheme
UASc=[1 1 1
    0 0 90/255
    212/255 102/255 0];
USAc=[1 1 1
    0 100/255 100/255
    212/255 102/255 0];
AUSc=[1 1 1
    0 0 90/255
    150/255 0 0];
Ec=[1 1 1
    0 0 0
    0 0 0];

% Stack the data
YUASR=[]; YUASL=[];
for i=1:37
    YUASR=[YUASR;[YSmR(i),r2realR(i)-YSmR(i),YUmR(i)-r2realR(i)]];
    YUASL=[YUASL;[YSmL(i),r2realL(i)-YSmL(i),YUmL(i)-r2realL(i)]];
end

% Monkey 1 right hemisphere
figure; hold on
a2=1:21;
HR=bar(a2,YUASR(a2,:),1,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(a2,YSmR(a2),YSseR(a2),'k.')
errorbar(a2,YUmR(a2),YUseR(a2),'k.')
title('Monkey 1 right hemisphere (Figure 2i)')

% Monkey 1 left hemisphere
figure; hold on
a2=1:21;
HL=bar(a2,YUASL(a2,:),1,'stack');
for k=1:3
    set(HL(k),'facecolor',UASc(k,:));
    set(HL(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(a2,YSmL(a2),YSseL(a2),'k.')
errorbar(a2,YUmL(a2),YUseL(a2),'k.')
title('Monkey 1 left hemisphere (supp)')

% Monkey 2 right hemisphere
figure; hold on
a2=22:37;
HR=bar(a2,YUASR(a2,:),1,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(a2,YSmR(a2),YSseR(a2),'k.')
errorbar(a2,YUmR(a2),YUseR(a2),'k.')
title('Monkey 2 right hemisphere (supp)')

% Monkey 2 left hemisphere
figure; hold on
a2=22:37;
HL=bar(a2,YUASL(a2,:),1,'stack');
for k=1:3
    set(HL(k),'facecolor',UASc(k,:));
    set(HL(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(a2,YSmL(a2),YSseL(a2),'k.')
errorbar(a2,YUmL(a2),YUseL(a2),'k.')
title('Monkey 2 left hemisphere (supp)')

% Mean bar figures
% Normalize so the upper bound is 1 for each rec. session
for i=1:37
    fctR=1/YUmR(i);
    YSmRsc(i)=YSmR(i)*fctR;
    r2realRsc(i)=r2realR(i)*fctR;
    YUmRsc(i)=YUmR(i)*fctR;
    fctL=1/YUmL(i);
    YSmLsc(i)=YSmL(i)*fctL;
    r2realLsc(i)=r2realL(i)*fctL;
    YUmLsc(i)=YUmL(i)*fctL;
end
% Stack the means
YM1Rm=[mean(YSmRsc(1:21)) mean(r2realRsc(1:21))-mean(YSmRsc(1:21)) mean(YUmRsc(1:21))-mean(r2realRsc(1:21))];
YM2Rm=[mean(YSmRsc(22:37)) mean(r2realRsc(22:37))-mean(YSmRsc(22:37)) mean(YUmRsc(22:37))-mean(r2realRsc(22:37))];
YM1Lm=[mean(YSmLsc(1:21)) mean(r2realLsc(1:21))-mean(YSmLsc(1:21)) mean(YUmLsc(1:21))-mean(r2realLsc(1:21))];
YM2Lm=[mean(YSmLsc(22:37)) mean(r2realLsc(22:37))-mean(YSmLsc(22:37)) mean(YUmLsc(22:37))-mean(r2realLsc(22:37))];

% Right hemisphere means
figure; hold on
HR=bar([1 2],[YM1Rm;YM2Rm],0.75,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
xlim([0 3])
title('Right hemisphere means, both monkeys')

% Left hemisphere means
figure; hold on
HR=bar([1 2],[YM1Lm;YM2Lm],0.75,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
xlim([0 3])
title('Left hemisphere means, both monkeys')