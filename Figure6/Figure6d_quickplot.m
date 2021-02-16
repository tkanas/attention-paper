% Plot figure 6d, kE vs area ratio

clear all
load areadata

% Plot the ratio of the computed area to the rectangle of dimensions
% 0.1x0.1 that the boundaries are computed in. That way units don't matter
% and it's less arbitrarily dependent on some of the parameters.
figure; plot(kEvec,areavec/0.01,'k','LineWidth',2);
xlabel('k_E')
ylabel('Area Ratio')
box off