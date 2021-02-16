% Compute and plot f-i curves for E and I populations

clear all

% units here are spikes per time constant tau (10 ms); multiplying by 100
% would give spikes per second
rEstart=0.0524; rEend=0.08;
rIstart=0.0573; rIend=0.14;

rEvec=[rEstart rEend];
rIvec=[rIstart rIend];

c=1; % overlap of noise to E and I populations, c=1 means they're 
% receiving the same noise

cei=0.05; % input correlation of pairs of neurons

sigE=0.3; sigI=0.35;
Vth=1; Vreset=0;

% Coupling
JE=1.5; JI=3;

% Transfer function
% erfcx is built-in MATLAB scaled complementary error function
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

% Units of time constant tau=1 (tau=1 is taken to be 10 ms in unit
% calculations)
taumax=10;
dt=0.005;
tau=-taumax:dt:taumax;

C11=zeros(length(rEvec),length(tau));

for j=1:length(rEvec)
    C11tmp=zeros(length(tau),1);
    rE=rEvec(j);
    rI=rIvec(j);
    
    % Solve for muE and muI from rE and rI
    % If first calculation, use muE,muI=1 as guesses, otherwise use
    % previously computed values as guesses
    if j==1
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,1,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,1,options);
    else
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,muE,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,muE,options);
    end
    
    IE0=muE+JE*rE-JI*rI;
    II0=muI+JE*rE-JI*rI;
   
    % LE=drE/dmueffE, LI=drI/dmueffI
    LE=rE^2*(sqrt(pi)/sigE)*(erfcx((-Vth+IE0)/sigE)-erfcx((-Vreset+IE0)/sigE));
    LI=rI^2*(sqrt(pi)/sigI)*(erfcx((-Vth+II0)/sigI)-erfcx((-Vreset+II0)/sigI));

    % drift matrix
    % first-order coefficients of linearization
    a_11=-1+JE*LE;
    a_22=-1-JI*LI;
    a_12=-JI*LE;
    a_21=JE*LI;
    A=-[a_11 a_12; a_21 a_22];
    [eV,eD]=eig(-A);

    % diffusion matrix
    D=[LE*sigE*sqrt(1-c) 0 LE*sigE*sqrt(c); 0 LI*sigI*sqrt(1-c) LI*sigI*sqrt(c)];

    % steady state variance matrix
    V=(det(A)*(D*D')+(A-trace(A)*eye(2))*(D*D')*(A-trace(A)*eye(2))')./(2*trace(A)*det(A));
    % correlation function
    for i=1:length(tau)
        if i<=length(tau)/2 % tau negative
            C=V*expm(A'*tau(i)); % Green's function
            C11tmp(i)=C(1);          
        else
            C=expm(-A*tau(i))*V; % Green's function
            C11tmp(i)=C(1);           
        end
    end
    C11(j,:)=cei*C11tmp;
end

% Plot with ms on x axis, (sp/s)^2 on y axis (multiplying rate by 100 gives
% sp/s, so multiplying by 10000 gives (sp/s)^2)
figure; plot(tau*10,C11*10000)
xlim([-25 25])
ylim([-0.5 3])
xlabel('Lag (ms)')
ylabel('Exc. autocov. function (sp/s^2)')
box off