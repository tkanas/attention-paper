clear all

sigE=0.3;
kE=1;
c=0.99;
cbar=0.01;
sigEbar=sqrt(cbar)*sigE;

JE=1.5; JI=3;
Vth=1;
Vreset=0;
rE0=0.01; rE1=0.1;

% Transfer function
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
options=optimset('Display','off','TolFun',1e-8);

muE0=fsolve(@(muE)fE(muE)-rE0,0,options);
muE1=fsolve(@(muE)fE(muE)-rE1,0,options);
dmu=muE1-muE0;

thetavec=0:0.01:0.5;
rEunvec=zeros(length(thetavec),1);
rEatvec=zeros(length(thetavec),1);
i=0;
for theta=thetavec
    i=i+1;
    % Unattended
    rEunvec(i)=fsolve(@(rE)fE(muE0+kE*theta)-rE,0,options);
    % Attended
    rEatvec(i)=fsolve(@(rE)fE(muE0+dmu+kE*theta)-rE,0,options);
end

figure; plot(thetavec,rEunvec*100,thetavec,rEatvec*100,'LineWidth',2)
xlim([0 0.3])
ylim([0 36])
xlabel('Stimulus intensity')
ylabel('Firing rate (sp/s)')
box off

theta=0.175;
x=0:0.001:0.5;
inputgauss=normpdf(x,0.175,sigEbar);

y=0:0.001:0.6;
IE0U=muE0+kE*theta;
rEU=fsolve(@(rE)fE(IE0U)-rE,0,options);
YU=rEU^2*(sqrt(pi)/sigE)*(exp((-Vth+IE0U)^2/(sigE^2))*erfc((-Vth+IE0U)/sigE)...
   -exp((-Vreset+IE0U)^2/(sigE^2))*erfc((-Vreset+IE0U)/sigE));
outputgaussU=normpdf(y,rEU,sigEbar*YU);

IE0A=muE0+dmu+kE*theta;
rEA=fsolve(@(rE)fE(IE0A)-rE,0,options);
YA=rEA^2*(sqrt(pi)/sigE)*(exp((-Vth+IE0A)^2/(sigE^2))*erfc((-Vth+IE0A)/sigE)...
   -exp((-Vreset+IE0A)^2/(sigE^2))*erfc((-Vreset+IE0A)/sigE));
outputgaussA=normpdf(y,rEA,sigEbar*YA);

figure; plot(x,inputgauss,'k','LineWidth',2);
yy=[-1 0 1];
hold on
plot(zeros(3,1)+theta-sigEbar,yy,zeros(3,1)+theta+sigEbar,yy)
hold off
xlim([0 0.3])
box off

figure; plot(outputgaussU,y*100,'LineWidth',2); hold on
xx=[-1 0 1];
plot(xx,(zeros(3,1)+rEU-sigEbar*YU)*100,xx,(zeros(3,1)+rEU+sigEbar*YU)*100);
plot(outputgaussA,y*100,'g','LineWidth',2)
plot(xx,(zeros(3,1)+rEA-sigEbar*YA)*100,xx,(zeros(3,1)+rEA+sigEbar*YA)*100)
ylim([0 36])
hold off
box off