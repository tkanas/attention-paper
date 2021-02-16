clear all

sigE=0.3;
sigI=0.35;

kE=1; kI=0.2;
c=0.99;
JE=1.5; JI=3;

Vth=1;
Vreset=0;

rEx0=0.0524; rEx1=0.08;
rIx0=0.0573; rIx1=0.14;

nmu=100;
drE=(rEx1-rEx0)/(nmu-1);
drI=(rIx1-rIx0)/(nmu-1);
rEvec=rEx0:drE:rEx1;
rIvec=rIx0:drI:rIx1;
nmu=length(rEvec);

% Transfer functions
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

mueffEvec=zeros(nmu,1);

fi=zeros(nmu,1);
VarEvec=zeros(nmu,1);
VarIvec=zeros(nmu,1);
GainEvec=zeros(nmu,1);
GainIvec=zeros(nmu,1);

i=0;
for rE=rEvec
    i=i+1;
    rI=rIvec(i);
    if i==1
        priorE=0;
        priorI=0;
    else
        priorE=muE;
        priorI=muI;
    end

    muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,priorE,options);
    muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,priorI,options);

    IE0=muE+JE*rE-JI*rI;
    II0=muI+JE*rE-JI*rI;
    mueffEvec(i)=IE0;

    LE=rE^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
    LI=rI^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

    GainE=(LE/(1+JI*LI-JE*LE))*(kE+JI*LI*(kE-kI));
    GainI=(LI/(1+JI*LI-JE*LE))*(kI+JE*LE*(kE-kI));

    M=-[-1+JE*LE -JI*LE;JE*LI -1-JI*LI];
    D=[LE*sigE*sqrt(1-c) 0 LE*sigE*sqrt(c); 0 LI*sigI*sqrt(1-c) LI*sigI*sqrt(c)];
    Q=inv(M)*D*D'*inv(M');
    CovEI=Q(1,2);
    VarE=Q(1,1);
    VarI=Q(2,2);

    VarEvec(i)=VarE;
    VarIvec(i)=VarI;
    GainEvec(i)=GainE;
    GainIvec(i)=GainI;
    fi(i)=(1/(VarE*VarI-CovEI^2))*(GainE^2*VarI+GainI^2*VarE-2*GainE*GainI*CovEI);
end

fiE=(GainEvec.^2)./VarEvec;
fiI=(GainIvec.^2)./VarIvec;

figure; plot(mueffEvec,fi,'k',mueffEvec,fiE,'k','LineWidth',2)
xlabel('Attention')
ylabel('Fisher Information')
box off