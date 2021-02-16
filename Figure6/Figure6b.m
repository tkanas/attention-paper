clear all
% compute gain along attentional path

kE=1;
kI=0;
K=[kE;kI];
JE=1.5;
JI=3;
theta=1;

Vth=1; Vreset=0;

sigE=0.3; sigI=0.35;

% Transfer functions
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

n=100;

rEstart=0.0524; rEend=0.08;
rIstart=0.0573; rIend=0.14;

rEvec=linspace(rEstart,rEend,n);
rIvec=linspace(rIstart,rIend,n);

gains=zeros(n,1);

i=0;
for rE=rEvec
    i=i+1;
    rI=rIvec(i);

    if i==1
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,1,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,1,options);
    else
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,muE,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,muI,options);
    end     
    
    IE0=muE+JE*rE-JI*rI;
    II0=muI+JE*rE-JI*rI;

    LE=rE^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
    LI=rI^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

    A=[-1+JE*LE,-JI*LI;JE*LE,-1-JI*LI];
    G=-inv(A)*K;
    gains(i)=LE*G(1);
end

Avec=linspace(0,1,n);
figure; hold on
plot(Avec,gains*100,'k','LineWidth',2)
plot(0,100*gains(1),'b.','MarkerSize',28);
plot(1,100*gains(end),'r.','MarkerSize',28);
xlabel('Attention')
ylabel('Gain')
xlim([0 1])
ylim([60 80])
box off
hold off