% Plot no-leave-out bar plots for all monkeys and hemispheres

clear all

% noleaveout_all contains pruned-down values to plot as stacked bar graphs. 
% It's more convenient than loading a bunch of unnecessary data.
% notation: S = shuffled, A = actual data, U = upper bound, 
% R = right hemisphere, L = left hemisphere, M1 = monkey 1 (datasets 1-21),
% M2 = monkey 2 (datasets 22-37)
load noleaveout_all

% Bar plot color scheme
% Face colors:
UASc=[1 1 1
    0 0 90/255
    212/255 102/255 0];
% Edge colors (black):
Ec=[1 1 1
    0 0 0
    0 0 0];

% Build up stacks
stackUASR=zeros(37,3); stackUASL=zeros(37,3);
for i=1:37
    stackUASR(i,:)=[YSmR(i),YAmR(i)-YSmR(i),YUmR(i)-YAmR(i)];
    stackUASL(i,:)=[YSmL(i),YAmL(i)-YSmL(i),YUmL(i)-YAmL(i)];
end

% Mean bar figs; normalize each bar so that the upper bound is 1
scaledSR=zeros(37,1);
scaledAR=zeros(37,1);
scaledUR=zeros(37,1);
scaledSL=zeros(37,1);
scaledAL=zeros(37,1);
scaledUL=zeros(37,1);
for i=1:37
    fctR=1/YUmR(i);
    scaledSR(i)=YSmR(i)*fctR;
    scaledAR(i)=YAmR(i)*fctR;
    scaledUR(i)=YUmR(i)*fctR;
    fctL=1/YUmL(i);
    scaledSL(i)=YSmL(i)*fctL;
    scaledAL(i)=YAmL(i)*fctL;
    scaledUL(i)=YUmL(i)*fctL;
end
stackedmeanM1R=[mean(scaledSR(1:21)) mean(scaledAR(1:21))-mean(scaledSR(1:21)) mean(scaledUR(1:21))-mean(scaledAR(1:21))];
stackedmeanM2R=[mean(scaledSR(22:37)) mean(scaledAR(22:37))-mean(scaledSR(22:37)) mean(scaledUR(22:37))-mean(scaledAR(22:37))];
stackedmeanM1L=[mean(scaledSL(1:21)) mean(scaledAL(1:21))-mean(scaledSL(1:21)) mean(scaledUL(1:21))-mean(scaledAL(1:21))];
stackedmeanM2L=[mean(scaledSL(22:37)) mean(scaledAL(22:37))-mean(scaledSL(22:37)) mean(scaledUL(22:37))-mean(scaledAL(22:37))];

% Plot the mean bar plots for monkeys 1 and 2 right hemisphere
figure; hold on
HR=bar([1 2],[stackedmeanM1R;stackedmeanM2R],0.75,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
xlim([0 3])
title('Averaged, scaled, RIGHT HEMISPHERE, monkey 1 (left), monkey 2 (right), Fig2h')

% Plot the mean bar plots for monkeys 1 and 2 left hemisphere
figure; hold on
HR=bar([1 2],[stackedmeanM1L;stackedmeanM2L],0.75,'stack');
for k=1:3
    set(HR(k),'facecolor',UASc(k,:));
    set(HR(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
xlim([0 3])
title('Averaged, scaled, LEFT HEMISPHERE, monkey 1 (left), monkey 2 (right), Fig2h')


% Individual bar plots
% Monkey 1 right hemisphere
figure; hold on
HR1=bar(1:21,stackUASR(1:21,:),1,'stack');
for k=1:3
    set(HR1(k),'facecolor',UASc(k,:));
    set(HR1(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(1:21,YSmR(1:21),YSseR(1:21),'k.')
errorbar(1:21,YUmR(1:21),YUseR(1:21),'k.')
title('Monkey 1 right hemisphere (Figure2g)')

% Monkey 2 right hemisphere
figure; hold on
HR2=bar(1:16,stackUASR(22:37,:),1,'stack');
for k=1:3
    set(HR2(k),'facecolor',UASc(k,:));
    set(HR2(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(1:16,YSmR(22:37),YSseR(22:37),'k.')
errorbar(1:16,YUmR(22:37),YUseR(22:37),'k.')
title('Monkey 2 right hemisphere (supp)')

% Monkey 1 left hemisphere
figure; hold on
HL1=bar(1:21,stackUASL(1:21,:),1,'stack');
for k=1:3
    set(HL1(k),'facecolor',UASc(k,:));
    set(HL1(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(1:21,YSmL(1:21),YSseL(1:21),'k.')
errorbar(1:21,YUmL(1:21),YUseL(1:21),'k.')
title('Monkey 1 left hemisphere (supp)')

% Monkey 2 left hemisphere
figure; hold on
HL2=bar(1:16,stackUASL(22:37,:),1,'stack');
for k=1:3
    set(HL2(k),'facecolor',UASc(k,:));
    set(HL2(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
errorbar(1:16,YSmL(22:37),YSseL(22:37),'k.')
errorbar(1:16,YUmL(22:37),YUseL(22:37),'k.')
title('Monkey 2 left hemisphere (supp)')