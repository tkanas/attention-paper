clear all
load covgainbdrydata

Covbdry=100*Covbdry;
Gainbdry=100*Gainbdry;
drIvec=drIvec*100;
rEend=8; rIend=14;

drIatt=rIend-100*rI0;
drEatt=rEend-100*rE0;

% Make a fit to get rid of artifact bumps (only arise because deltacov=0
% was computed using minima)
Fitcov = polyfit(drIvec',Covbdry,2);
Fitgain = polyfit(drIvec',Gainbdry,2);

figure; hold on
plot(drIvec,polyval(Fitcov,drIvec),'k','LineWidth',2)
plot(drIvec,polyval(Fitgain,drIvec),'k','LineWidth',2);

l=line([0 drIatt],[0 drEatt],'LineWidth',2);
plot(0,0,'b.','MarkerSize',28);
plot(drIatt,drEatt,'r.','MarkerSize',28);
plot(0:10,0:10,'k--')
hold off
xlim([0 10])
ylim([0 10])
xlabel('drI')
ylabel('drE')
box off