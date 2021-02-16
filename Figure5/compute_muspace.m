% Compute the boundary that separates delta-mu-space into delta-cov<0 and
% delta-cov>0 region

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

% Unattended rates
rE0=0.0524;
rI0=0.0573;

% Unattended muE, muI, mueffE, mueffI
muE0=fsolve(@(muE)fE(JE*rE0-JI*rI0+muE)-rE0,1,options);
muI0=fsolve(@(muI)fI(JE*rE0-JI*rI0+muI)-rI0,1,options);
IE0=muE0+JE*rE0-JI*rI0;
II0=muI0+JE*rE0-JI*rI0;

% set up for iterating through muE-muI space
dmuE0=0; dmuI0=0;
dmuEend=0.4; dmuIend=0.4;
n=100;
dmuEvec=linspace(dmuE0,dmuEend,n);
dmuIvec=linspace(dmuI0,dmuIend,n);

% LE=drE/dmueffE, LI=drI/dmueffI
LE0=rE0^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
LI0=rI0^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

% Unattended covariances
CovE0=((LE0*(JI*LI0*(sigE-sigI)+sigE))/(1+JI*LI0-JE*LE0))^2;
CovI0=((LI0*(JE*LE0*(sigE-sigI)+sigI))/(1+JI*LI0-JE*LE0))^2;

CovEEbdry=zeros(n,1); % boundary between dcovEE<0 and dcovEE>0
dcovEE=zeros(n);

i=0;
for dmuI=dmuIvec
    i=i+1
    muI1=muI0+dmuI;
    j=0;
    for dmuE=dmuEvec
        j=j+1;
        muE1=muE0+dmuE;
        if i==1 && j==1
            r=fsolve(@(r)[fE(JE*r(1)-JI*r(2)+muE1)-r(1);fI(JE*r(1)-JI*r(2)+muI1)-r(2)],[0;0],options);
            rE1=r(1); rI1=r(2);
        else
            r=fsolve(@(r)[fE(JE*r(1)-JI*r(2)+muE1)-r(1);fI(JE*r(1)-JI*r(2)+muI1)-r(2)],[rE1;rI1],options);
            rE1=r(1); rI1=r(2);
        end
        
        IE1=muE1+JE*rE1-JI*rI1;
        II1=muI1+JE*rE1-JI*rI1;
         
        LE1=rE1^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE1)/sigE)-erfcx((-Vreset+IE1)/sigE));
        LI1=rI1^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II1)/sigI)-erfcx((-Vreset+II1)/sigI));
        
        lambda=-1-JI*LI1+JE*LE1;
        
        CovE1=((LE1*(JI*LI1*(sigE-sigI)+sigE))/(1+JI*LI1-JE*LE1))^2;
        CovI1=((LI1*(JE*LE1*(sigE-sigI)+sigI))/(1+JI*LI1-JE*LE1))^2;
        
        if lambda<0 % stable
            dcovEE(j,i)=CovE1-CovE0;
        else
            dcovEE(j,i)=NaN;
        end
    end
    diffveccovE=abs(dcovEE(:,i));
    mindifcovE=min(diffveccovE);
    bestrind=find(diffveccovE==mindifcovE);
    if ~isnan(mindifcovE)
        CovEEbdry(i)=dmuEvec(bestrind);
    else
        CovEEbdry(i)=NaN;
    end
end

% Attended mu
rE1=0.08;
rI1=0.14;
muE1=fsolve(@(muE)fE(JE*rE1-JI*rI1+muE)-rE1,1,options);
muI1=fsolve(@(muI)fI(JE*rE1-JI*rI1+muI)-rI1,1,options);
dmuIatt=muI1-muI0;
dmuEatt=muE1-muE0;

save dmudata1