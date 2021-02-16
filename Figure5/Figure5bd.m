% compute population variance and eigenvalues along attentional path

clear all

sigE=0.3; sigI=0.35;
JE=1.5; JI=3;
Vth=1;
Vreset=0;

% Transfer functions
% erfcx is built-in MATLAB scaled complementary error function
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

% Window of rE-rI space to compute VarE over
rEstart=0.0524; rEend=0.08;
rIstart=0.0573; rIend=0.14;

nmu=50;
rEvec=linspace(rEstart,rEend,nmu);
rIvec=linspace(rIstart,rIend,nmu);

VE=zeros(nmu,1);
eigs=zeros(nmu,1);
mueffEvec=zeros(nmu,1);
mueffIvec=zeros(nmu,1);

i=0;
for rE=rEvec
    i=i+1;
    rE
    rI=rIvec(i);

    % Solve for muE and muI from rE and rI
    % If first calculation, use muE,muI=1 as guesses, otherwise use
    % previously computed values as guesses
    if i==1
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,1,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,1,options);
    else
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,muE,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,muE,options);
    end

    % mu-effective
    IE0=muE+JE*rE-JI*rI;
    II0=muI+JE*rE-JI*rI;
    mueffEvec(i)=IE0;
    mueffIvec(i)=II0;

    % LE=drE/dmueffE, LI=drI/dmueffI
    LE=rE^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
    LI=rI^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

    VE(i)=((LE*(JI*LI*(sigE-sigI)+sigE))/(1+JI*LI-JE*LE))^2;
    eigs(i)=-1-JI*LI+JE*LE;
end

attpath=linspace(0,1,nmu);
figure; plot(attpath,VE*10000,'k')
hold on
plot(0,VE(1)*10000,'g.','MarkerSize',28)
plot(1,VE(end)*10000,'m.','MarkerSize',28)
hold off
ylabel('E population variance')
xlabel('Attention')
ylim([20 50])
box off

figure; plot(attpath,eigs,'k')
hold on
plot(0,eigs(1),'g.','MarkerSize',28)
plot(1,eigs(end),'m.','MarkerSize',28)
hold off
ylabel('Eigenvalue')
box off
xlabel('Attention')