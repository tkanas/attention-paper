clear all
load surrconvpf

figure; hold on
for i=1:ndata
    plot(squeeze(extrvec(1,i,:)),squeeze(extrvec(2,i,:)),'k')
end
plot(perfvec(1,:),perfvec(2,:),'r.','MarkerSize',10)
hold off
xlabel('Number of Samples')
ylabel('r^2')