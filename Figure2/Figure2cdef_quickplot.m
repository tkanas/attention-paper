% Plot scatter plots for Fig. 2 b,c,d
% Data for recording session 10, right hemisphere, no-leave-out

clear all
load scatterplot_data

% comprealR contains predicted (first column) and actual (second column)
% covariance values for unaltered dataset
figure; plot(comprealR(:,2),comprealR(:,1),'k.')
hold on
xlabel('Actual Covariance')
ylabel('Predicted Covariance')
title('Data')
xlim([-2 6])
ylim([-2 6])
plot(-2:6,-2:6,'k')
box off
hold off

% compshufR contains predicted (first column) and actual (second column)
% covariance values for shuffled dataset
figure; plot(compshufR(:,2),compshufR(:,1),'k.')
hold on
xlabel('Actual Covariance')
ylabel('Predicted Covariance')
title('Shuffled')
xlim([-2 6])
ylim([-2 6])
plot(-2:6,-2:6,'k')
box off
hold off

% compsurrR contains predicted (first column) and actual (second column)
% covariance values for upper bound dataset
figure; plot(compsurrR(:,2),compsurrR(:,1),'k.')
hold on
xlabel('Actual Covariance')
ylabel('Predicted Covariance')
title('Upper bound')
xlim([-2 6])
ylim([-2 6])
plot(-2:6,-2:6,'k')
box off
hold off

% color scheme for bar plot
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

% compute r^2 for this example
r2real=corrcoef(comprealR);
r2real=r2real(1,2);
r2shuf=corrcoef(compshufR);
r2shuf=r2shuf(1,2);
r2ub=corrcoef(compsurrR);
r2ub=r2ub(1,2);

% bar values to be plotted
Y=[r2shuf,r2real-r2shuf,r2ub-r2real];
Y(2,3)=0; % Silly matlab trick: have to pad with zeros to plot single column

figure; hold on
H=bar([1 2],Y,0.5,'stack');
for k=1:3
    set(H(k),'facecolor',UASc(k,:));
    set(H(k),'edgecolor',Ec(k,:));
end
set(gca,'FontSize',14);
xlim([0 2])