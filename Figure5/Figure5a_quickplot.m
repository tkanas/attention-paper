% Plot attentional path through contours of VarE

%clear all

load varvalues
cei=0.05;
VE=VE*cei;

figure; [C,h]=contour(rIvec*100,rEvec*100,VE*10000);
h.LevelList=[10 20 30 40 50 60 70 80];
h.LevelList=h.LevelList*cei;
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