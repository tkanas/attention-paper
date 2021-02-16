clear all
load gains2D

figure; [C,h]=contour(100*rIvec,100*rEvec,100*gains);
h.LevelList=40:10:110;
hold on
l=line([5.73 14],[5.24 8],'LineWidth',2);
plot(5.73,5.24,'b.','MarkerSize',28);
plot(14,8,'r.','MarkerSize',28);
xlabel('rI')
ylabel('rE')
box off
colormap winter
colorbar
hold off