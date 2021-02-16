clear all


JE=1.5; JI=3;
rE0=0.0524; rI0=0.0573;
drE0=0; drI0=0;
drEend=0.1; drIend=0.1;

n=100;
drEvec=linspace(drE0,drEend,n);
drIvec=linspace(drI0,drIend,n);

Vth=1; Vreset=0;
kE=1; kI=0;
c=1;
sigE=0.3; sigI=0.35;

% Transfer functions
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

muE0=fsolve(@(muE)fE(JE*rE0-JI*rI0+muE)-rE0,0,options);
muI0=fsolve(@(muI)fI(JE*rE0-JI*rI0+muI)-rI0,0,options);
IE0=JE*rE0-JI*rI0+muE0;
II0=JE*rE0-JI*rI0+muI0;
LE0=rE0^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
LI0=rI0^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

Cov0=((LE0*(JI*LI0*(sigE-sigI)+sigE))/(1+JI*LI0-JE*LE0))^2;
G0=(LE0*(kE+JI*LI0*(kE-kI)))/(1+JI*LI0-JE*LE0);

Covbdry=zeros(n,1);
Gainbdry=zeros(n,1);

differenzcov=zeros(n);
differenzgain=zeros(n);
gainvec=zeros(n);
covvec=zeros(n);

i=0;
for drI=drIvec
    i=i+1
    rI1=rI0+drI;
    j=0;
    for drE=drEvec
        j=j+1;
        rE1=rE0+drE;
        if i==1 && j==1
            [muE1,~,efE]=fsolve(@(muE)fE(JE*rE1-JI*rI1+muE)-rE1,0,options);
            [muI1,~,efI]=fsolve(@(muI)fI(JE*rE1-JI*rI1+muI)-rI1,0,options);
            foqi=[muE1,muI1];
        elseif j>1 && j<=length(drEvec)
            [muE1,~,efE]=fsolve(@(muE)fE(JE*rE1-JI*rI1+muE)-rE1,muE1,options);
            [muI1,~,efI]=fsolve(@(muI)fI(JE*rE1-JI*rI1+muI)-rI1,muI1,options);
        elseif i>1 && j==1
            [muE1,~,efE]=fsolve(@(muE)fE(JE*rE1-JI*rI1+muE)-rE1,foqi(1),options);
            [muI1,~,efI]=fsolve(@(muI)fI(JE*rE1-JI*rI1+muI)-rI1,foqi(2),options);
            foqi=[muE1,muI1];
        else
            error('missed condition');
        end
        if efE~=1 || efI~=1
            error('not solved')
        end
        IE1=muE1+JE*rE1-JI*rI1;
        II1=muI1+JE*rE1-JI*rI1;
        
        LE1=rE1^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE1)/sigE)-erfcx((-Vreset+IE1)/sigE));
        LI1=rI1^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II1)/sigI)-erfcx((-Vreset+II1)/sigI));
        
        Cov1=LE1^2*(JI*LI1*(sigE-sigI)+sigE)^2/(1+JI*LI1-JE*LE1)^2;
        GE1=(LE1*(kE+JI*LI1*(kE-kI)))/(1+JI*LI1-JE*LE1);
        if -1-JI*LI1+JE*LE1<0 % stable
            differenzcov(j,i)=Cov1-Cov0;
            differenzgain(j,i)=GE1-G0;
            gainvec(j,i)=GE1;
            covvec(j,i)=Cov1;
        else
            differenzcov(j,i)=NaN;
            differenzgain(j,i)=NaN;
            gainvec(j,i)=NaN;
            covvec(j,i)=NaN;
        end
    end
    diffveccov=abs(differenzcov(:,i));
    mindifcov=min(diffveccov);
    bestrind=find(diffveccov==mindifcov);
    if ~isnan(mindifcov)
        Covbdry(i)=drEvec(bestrind);
    else
        Covbdry(i)=NaN;
    end
    diffvecgain=abs(differenzgain(:,i));
    mindifgain=min(diffvecgain);
    bestrind=find(diffvecgain==mindifgain);
    if ~isnan(mindifgain)
        Gainbdry(i)=drEvec(bestrind);
    else
        Gainbdry(i)=NaN;
    end
end


save covgainbdrydata1