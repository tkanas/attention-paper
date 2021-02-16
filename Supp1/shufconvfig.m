clear all
load shufconvdata
load noleaveout_all

figure; plot(Nvec,r2vec,'k')
hold on
plot(nunits,YSmR,'r.','MarkerSize',20)

set(gca,'FontSize',14);
xlabel('Number of units')
ylabel('r^2')
box off