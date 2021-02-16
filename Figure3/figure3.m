% Create figure 3: cov plot
clear all

figure;

load 2Dcovnew
subplot(2,2,1)
imagesc(100*rIvec,100*rEvec,log10(10000*CREREc))
h=colorbar;
set(h,'YScale','log')
shading interp
set(gca,'FontSize',14)
set(gca,'YDir','normal')
set(gca,'CLim',[0 2])
colormap(gray)
hold on
l=line([5.73 14],[5.24 8],'LineWidth',2);
plot(5.73,5.24,'b.','LineWidth',2);
plot(14,8,'r.','LineWidth',2);
%plot(100*rIvec,100*rEbdry,'r','LineWidth',2)
hold off
xlabel('rI')
ylabel('rE')
xlim([1 18]);
ylim([1 18]);

load 1dcovplot
subplot(2,2,2)
x=linspace(0,1,100);
plot(x,CREREc*10000,'k','LineWidth',2)
%ylim([0.71 1.1])
set(gca,'FontSize',14)
xlabel('A')
ylabel('Pop Variance (spikes/s^2)')

subplot(2,2,3)
plot(x,-null,'k','LineWidth',2)
ylim([-2.2 -1.5])
set(gca,'FontSize',14)
xlabel('A')
ylabel('Eigenvalue')

load autocovun
subplot(2,2,4)
plot(100*tau,10000*C11,'LineWidth',2)
hold on
load autocovat
plot(100*tau,10000*C11,'r','LineWidth',2)
hold off
set(gca,'FontSize',14)
xlim([-250 250])
%legend('Unattended','Attended')
xlabel('Lag (ms)')
ylabel('Exc. autocov. function (spikes/s^2)')