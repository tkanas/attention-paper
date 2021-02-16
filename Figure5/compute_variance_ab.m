% compute population variance in rE-rI space

clear all

sigE=0.3; sigI=0.35;

alpha = .9; beta = 0.9;
%JE=1.5; JI=3;
JEE=1.5; JIE=alpha*JEE;
JII = 3; JEI = beta*JII;
Vth=1;
Vreset=0;

% Transfer functions
erfcxfxn=@(z)erfcx(z);
fE = @(x) 1/(sqrt(pi)*integral(erfcxfxn,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(erfcxfxn,(-Vth+x)/sigI,(-Vreset+x)/sigI));
opts=optimset('Display','off');

% Window of rE-rI space to compute VarE over
rEx0=0.001; rIx0=0.001;
rEx1=0.18; rIx1=0.18;

nmu=100;
drE=(rEx1-rEx0)/(nmu-1);
drI=(rIx1-rIx0)/(nmu-1);
rEvec=rEx0:drE:rEx1;
rIvec=rIx0:drI:rIx1;

%rEvec=[0.0524,0.08]; rIvec=[0.0573,0.14];

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
        [muE0,f,efE]=fsolve(@(muE)fE(JEE*rE-JEI*rI+muE)-rE,10,opts);
        [muI0,f,efI]=fsolve(@(muI)fI(JIE*rE-JII*rI+muI)-rI,10,opts);
        if efE~=1 || efI~=1
            error('not solved')
        end
        
        % mu-effective
        IE0=JEE*rE-JEI*rI+muE0;
        II0=JIE*rE-JII*rI+muI0;
        mueffEvec(i,j)=IE0;
        mueffIvec(i,j)=II0;

        % LE=drE/dmueffE, LI=drI/dmueffI
        Y0=rE^2*(sqrt(pi)/sigE)*(exp((-Vth+IE0)^2/(sigE^2))*erfc((-Vth+IE0)/sigE)...
            -exp((-Vreset+IE0)^2/(sigE^2))*erfc((-Vreset+IE0)/sigE));
        Z0=rI^2*(sqrt(pi)/sigI)*(exp((-Vth+II0)^2/(sigI^2))*erfc((-Vth+II0)/sigI)...
            -exp((-Vreset+II0)^2/(sigI^2))*erfc((-Vreset+II0)/sigI));

        VE(i,j)=Y0^2*(JEI*Z0*sigI-JII*Z0*sigE-sigE)^2/(JEE*JII*Y0*Z0-JEI*JIE*Y0*Z0+JEE*Y0-JII*Z0-1)^2;

    end
end

save varvalues1