clear all
close all

sigE=0.3;
sigI=0.35;
kE=1;
kI=0.2;
c=0.99;
cbar=0.01;
sigEbar=sqrt(cbar)*sigE;
sigIbar=sqrt(cbar)*sigI;
theta0=0.175;

JE=1.5; JI=3;
Vth=1;
Vreset=0;
%rE0=0.01; rE1=0.1;
%rI0=0.01; rI1=0.1;

rE0=0.0524; rE1=0.08;
rI0=0.0573; rI1=0.14;

experf = @(z) exp(z.^2).*erfc(z);
fE = @(x) 1/(sqrt(pi)*integral(experf,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(experf,(-Vth+x)/sigI,(-Vreset+x)/sigI));
opts=optimset('Display','off');

muE0=fsolve(@(muE)fE(muE+kE*theta0+JE*rE0-JI*rI0)-rE0,0,opts);
muE1=fsolve(@(muE)fE(muE+kE*theta0+JE*rE1-JI*rI1)-rE1,0,opts);
dmuE=muE1-muE0;
muI0=fsolve(@(muI)fI(muI+kI*theta0+JE*rE0-JI*rI0)-rI0,0,opts);
muI1=fsolve(@(muI)fI(muI+kI*theta0+JE*rE1-JI*rI1)-rI1,0,opts);
dmuI=muI1-muI0;


thetavec=0.1:0.005:0.25;
rEunvec=zeros(length(thetavec),1);
rEatvec=zeros(length(thetavec),1);
i=0;
for theta=thetavec
    i=i+1;
    rU=fsolve(@(r)[fE(muE0+kE*theta+JE*r(1)-JI*r(2))-r(1);fI(muI0+kI*theta+JE*r(1)-JI*r(2))-r(2)],[0;0],opts);
    rEunvec(i)=rU(1);
    rA=fsolve(@(r)[fE(muE0+dmuE+kE*theta+JE*r(1)-JI*r(2))-r(1);fI(muI0+dmuI+kI*theta+JE*r(1)-JI*r(2))-r(2)],[0;0],opts);
    rEatvec(i)=rA(1);
end

figure; plot(thetavec,rEunvec*100,thetavec,rEatvec*100,'LineWidth',2)
xlim([0.1 0.25])
ylim([2 15])
box off
xlabel('Stimulus Intensity')
ylabel('Firing rate (sp/s)')

theta=theta0;

y=0:0.001:0.2;
rU=[rE0;rI0];
IE0U=muE0+kE*theta+JE*rU(1)-JI*rU(2);
II0U=muI0+kI*theta+JE*rU(1)-JI*rU(2);
rEU=rU(1);
rIU=rU(2);
YU=rEU^2*(sqrt(pi)/sigE)*(exp((-Vth+IE0U)^2/(sigE^2))*erfc((-Vth+IE0U)/sigE)...
   -exp((-Vreset+IE0U)^2/(sigE^2))*erfc((-Vreset+IE0U)/sigE));
ZU=rIU^2*(sqrt(pi)/sigI)*(exp((-Vth+II0U)^2/(sigI^2))*erfc((-Vth+II0U)/sigI)...
   -exp((-Vreset+II0U)^2/(sigI^2))*erfc((-Vreset+II0U)/sigI));

rA=[rE1;rI1];
IE0A=muE0+dmuE+kE*theta+JE*rA(1)-JI*rA(2);
II0A=muI0+dmuI+kI*theta+JE*rA(1)-JI*rA(2);
rEA=rA(1);
rIA=rA(2);
YA=rEA^2*(sqrt(pi)/sigE)*(exp((-Vth+IE0A)^2/(sigE^2))*erfc((-Vth+IE0A)/sigE)...
   -exp((-Vreset+IE0A)^2/(sigE^2))*erfc((-Vreset+IE0A)/sigE));
ZA=rIA^2*(sqrt(pi)/sigI)*(exp((-Vth+II0A)^2/(sigI^2))*erfc((-Vth+II0A)/sigI)...
   -exp((-Vreset+II0A)^2/(sigI^2))*erfc((-Vreset+II0A)/sigI));

%VarEU=(YU^2/(1+JI*ZU-JE*YU)^2)*(JI*ZU*(sigE-sigI)+sigE)^2;
%VarEA=(YA^2/(1+JI*ZA-JE*YA)^2)*(JI*ZA*(sigE-sigI)+sigE)^2;

%VarIU=(ZU^2/(1+JI*ZU-JE*YU)^2)*(JE*YU*(sigE-sigI)+sigI)^2;

tt=0.01;
MU=-[-1+JE*YU -JI*YU YU*sigE*sqrt(1-c)*sqrt(cbar) 0 YU*sigE*sqrt(c)*sqrt(cbar);
    JE*ZU -1-JI*ZU 0 ZU*sigI*sqrt(1-c)*sqrt(cbar) ZU*sigI*sqrt(c)*sqrt(cbar);
    0 0 -1/tt 0 0;
    0 0 0 -1/tt 0;
    0 0 0 0 -1/tt];
MA=-[-1+JE*YA -JI*YA YA*sigE*sqrt(1-c)*sqrt(cbar) 0 YA*sigE*sqrt(c)*sqrt(cbar);
    JE*ZA -1-JI*ZA 0 ZA*sigI*sqrt(1-c)*sqrt(cbar) ZA*sigI*sqrt(c)*sqrt(cbar);
    0 0 -1/tt 0 0;
    0 0 0 -1/tt 0;
    0 0 0 0 -1/tt];
D=[0 0 0;
    0 0 0;
    1/tt 0 0;
    0 1/tt 0;
    0 0 1/tt];
SU=inv(MU)*D*D'*inv(MU');
SA=inv(MA)*D*D'*inv(MA');

VarEU=SU(1,1);
VarEA=SA(1,1);

outputgaussU=normpdf(y,rEU,sqrt(VarEU));
outputgaussA=normpdf(y,rEA,sqrt(VarEA));

x=0.1:0.001:0.25;
inputvarU=JE^2*SU(1,1)-2*JE*JI*SU(1,2)+JI^2*SU(2,2)+2*sigE*sqrt(cbar)...
    *(JE*sqrt(1-c)*SU(1,3)+JE*sqrt(c)*SU(1,5)-JI*sqrt(1-c)*SU(2,3)-JI*sqrt(c)*SU(2,5))...
    +sigE^2*cbar*(1-c)*SU(3,3)+sigE^2*cbar*c*SU(5,5);
inputvarA=JE^2*SA(1,1)-2*JE*JI*SA(1,2)+JI^2*SA(2,2)+2*sigE*sqrt(cbar)...
    *(JE*sqrt(1-c)*SA(1,3)+JE*sqrt(c)*SA(1,5)-JI*sqrt(1-c)*SA(2,3)-JI*sqrt(c)*SA(2,5))...
    +sigE^2*cbar*(1-c)*SA(3,3)+sigE^2*cbar*c*SA(5,5);
inputgaussU=normpdf(x,theta,sqrt(inputvarU));
inputgaussA=normpdf(x,theta,sqrt(inputvarA));

figure; plot(outputgaussU,y*100,'LineWidth',2); hold on
xx=[-1 0 1];
plot(xx,(zeros(3,1)+rEU-sqrt(VarEU))*100,xx,(zeros(3,1)+rEU+sqrt(VarEU))*100);
plot(outputgaussA,y*100,'g','LineWidth',2)
plot(xx,(zeros(3,1)+rEA-sqrt(VarEA))*100,xx,(zeros(3,1)+rEA+sqrt(VarEA))*100)
hold off
ylim([2 15])
box off

figure; plot(x,inputgaussU);
yy=xx;
hold on
plot(zeros(3,1)+theta-sqrt(inputvarU),yy,zeros(3,1)+theta+sqrt(inputvarU),yy)
plot(x,inputgaussA);
plot(zeros(3,1)+theta-sqrt(inputvarA),yy,zeros(3,1)+theta+sqrt(inputvarA),yy)
hold off
xlim([0.1 0.25])
box off