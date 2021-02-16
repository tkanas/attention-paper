clear all
load dmudata

% Make a fit to get rid of artifact bumps (only arise because deltacov=0
% was computed using minima)
Fit = polyfit(dmuIvec',CovEEbdry,2);

figure; hold on
plot(dmuIvec,polyval(Fit,dmuIvec),'k','LineWidth',2)
l=line([0 dmuIatt],[0 dmuEatt],'LineWidth',2);
plot(0,0,'b.','MarkerSize',28);
plot(dmuIatt,dmuEatt,'r.','MarkerSize',28);
plot(0:0.01:0.4,0:0.01:0.4,'k--')
hold off
xlim([0 0.4])
ylim([0 0.4])
xlabel('dmuI')
ylabel('dmuE')
box off