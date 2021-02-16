clear all

Vth=1;
Vreset=0;
sigE=0.3;
sigI=0.35;

JEE=1.5;
JEI=1.5;
JIE=2.25;
JII=3;

erfcxfxn=@(z)erfcx(z);
fE = @(x) 1/(sqrt(pi)*integral(erfcxfxn,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(erfcxfxn,(-Vth+x)/sigI,(-Vreset+x)/sigI));
opts=optimset('Display','off');

rE=0.08;
rI=0.14;

[muI0,f,efI]=fsolve(@(muI)fI(JIE*rE-JII*rI+muI)-rI,10,opts);