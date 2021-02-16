% compute population variance in rE-rI space

clear all

sigE=0.3; sigI=0.35;
JE=1.5; JI=3;
Vth=1;
Vreset=0;

% Transfer functions
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

% Window of rE-rI space to compute VarE over
rEx0=0.001; rIx0=0.001;
rEx1=0.18; rIx1=0.18;

nmu=100;
drE=(rEx1-rEx0)/(nmu-1);
drI=(rIx1-rIx0)/(nmu-1);
rEvec=rEx0:drE:rEx1;
rIvec=rIx0:drI:rIx1;

nmu=length(rEvec);

VE=zeros(nmu);
mueffEvec=zeros(nmu);
mueffIvec=zeros(nmu);

i=0;
for rE=rEvec
    i=i+1;
    rE
    j=0;
    for rI=rIvec
        j=j+1;

        % Solve for muE and muI from rE and rI
        % If first calculation, use muE,muI=1 as guesses, otherwise use
        % previously computed values as guesses
        if i==1 && j==1
            muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,1,options);
            muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,1,options);
        else
            muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,muE,options);
            muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,muI,options);
        end

        % mu-effective
        IE0=muE+JE*rE-JI*rI;
        II0=muI+JE*rE-JI*rI;
        mueffEvec(i,j)=IE0;
        mueffIvec(i,j)=II0;

        % LE=drE/dmueffE, LI=drI/dmueffI
        LE=rE^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
        LI=rI^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

        VE(i,j)=((LE*(JI*LI*(sigE-sigI)+sigE))/(1+JI*LI-JE*LE))^2;

    end
end

save varvalues1