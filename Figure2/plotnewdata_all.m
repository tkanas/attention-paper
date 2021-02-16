clear all
load sname

grightall = [];
gleftall = [];
raterightall_ratio = [];
rateleftall_ratio = [];
num_neur_right=zeros(37,1);
num_neur_left = zeros(37,1);
for ndi=1:37
    n_neur = sum(~isnan(xvecallrealdataR(ndi,:)));
    grightall = [grightall xvecallrealdataR(ndi,1:n_neur)];
    num_neur_right(ndi) = n_neur;
    raterightall_ratio = [raterightall_ratio allratesatR(ndi,1:n_neur)./allratesunR(ndi,1:n_neur)];

    n_neur = sum(~isnan(xvecallrealdataL(ndi,:)));
    gleftall = [gleftall xvecallrealdataL(ndi,1:n_neur)];
    num_neur_left(ndi) = n_neur;
    rateleftall_ratio = [rateleftall_ratio allratesatL(ndi,1:n_neur)./allratesunL(ndi,1:n_neur)];
end

figure; hold on
[n1, xout1] = hist(grightall,30);
bar(xout1,n1,'r');
[n2, xout2] = hist(gleftall,30);
bar(xout2,n2,'g');
legend('Right hemisphere, n=1979','Left hemisphere, n=1181')
xlabel('Value of g')
ylabel('Occurrences')
hold off

%{
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
%}

figure; hold on
plot(raterightall_ratio, grightall,'ro')
plot(rateleftall_ratio, gleftall,'go')
hold off
xlabel('rA/rU')
ylabel('g')
title('ratrate')
