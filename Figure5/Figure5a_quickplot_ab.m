% Plot attentional path through contours of VarE

clear all

load varvaluesa09b09
cei=0.05;
VE=VE*cei;

figure; [C,h]=contour(rIvec*100,rEvec*100,VE*10000);
h.LevelList=[5 10 15 20 30 40 50 60 70 80];
% For alpha = 1.5, beta = 1.1:
%h.LevelList = linspace(1,20,8);
% For alpha = 1.5, beta = 0.8:
%h.LevelList = linspace(30,100,8);
% Fpr alpha = 0.5, beta = 1.2:
%h.LevelList = linspace(1,60,8);
% For alpha, beta = 0.9:
h.LevelList = linspace(30,130,8);
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