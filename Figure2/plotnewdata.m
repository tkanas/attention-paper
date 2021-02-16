clear all
load scatterplot_data

gright=xvecallrealdataR(10,1:39);
gleft=xvecallrealdataL(10,1:14);

figure; hold on
[n1, xout1] = hist(gright,10);
bar(xout1,n1,'r');
[n2, xout2] = hist(gleft,10);
bar(xout2,n2,'g');
legend('Right hemisphere, n=39','Left hemisphere, n=14')
xlabel('Value of g')
ylabel('Occurrences')
hold off

rrightunatt = allratesunR(10,1:39);
rrightatt = allratesatR(10,1:39);
rleftunatt = allratesunL(10,1:14);
rleftatt = allratesatL(10,1:14);

deltarright = rrightatt - rrightunatt;
deltarleft = rleftatt - rleftunatt;

ratrright = rrightatt./rrightunatt;
rratleft = rleftatt./rleftunatt;

figure; hold on
plot(deltarright,gright,'ro')
plot(deltarleft,gleft,'go')
hold off
xlabel('rA-rU')
ylabel('g')
title('deltarate')

figure; hold on
plot(ratrright,gright,'ro')
plot(rratleft,gleft,'go')
hold off
xlabel('rA/rU')
ylabel('g')
title('ratrate')

