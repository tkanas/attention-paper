% Compute and plot f-i curves for E and I populations

clear all

% Start and endpoints of attentional path in firing rate space
rEstart=0.0524; rEend=0.08;
rIstart=0.0573; rIend=0.14;

% Attentional path in firing rate space
rEvec=linspace(rEstart,rEend,30);
rIvec=linspace(rIstart,rIend,30);
drE=rEvec(2)-rEvec(1);
drI=rIvec(2)-rIvec(1);

% Extend FI curve around attentional path
% extra points below start of attentional path for E,I
nextralowE=floor(rEstart/drE);
nextralowI=floor(rIstart/drI);
% extra points above end of attentional path for E,I
nextrahighE=ceil((0.2-rEend)/drE);
nextrahighI=ceil((0.25-rIend)/drI);
% extended intervals for FI curves
rEextendedvec=rEstart-drE*nextralowE:drE:rEend+drE*nextrahighE;
rIextendedvec=rIstart-drI*nextralowI:drI:rIend+drI*nextrahighI;
nrE=length(rEextendedvec);
nrI=length(rIextendedvec);

sigE=0.3; sigI=0.35;
Vth=1; Vreset=0;

JE=1.5; JI=3; 

% Transfer function
fE = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigE,(-Vreset+x)/sigE)); 
fI = @(x) 1/(sqrt(pi)*integral(@erfcx,(-Vth+x)/sigI,(-Vreset+x)/sigI));
options=optimset('Display','off','TolFun',1e-8);

% extended mu-effectives
mueffEvec_ext=zeros(nrE,1);
mueffIvec_ext=zeros(nrI,1);

% mu-effectives along attentional path
mueffEvec=zeros(length(rEvec),1);
mueffIvec=zeros(length(rIvec),1);

% First compute attentional path
for j=1:length(rEvec)
    rE=rEvec(j);
    rI=rIvec(j);

    if j==1
        prior=-1;
    else
        prior(1)=muE;
        prior(2)=muI;
    end
    %[muE,muI]=mufromr(prior,1,rE,rI,JEE,JEI,JIE,JII,sigE,sigI,tauE,tauI,Vth,Vreset,fE,fI);
    if j==1
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,1,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,1,options);
    else
        muE=fsolve(@(muE)fE(JE*rE-JI*rI+muE)-rE,muE,options);
        muI=fsolve(@(muI)fI(JE*rE-JI*rI+muI)-rI,muI,options);
    end
    
    IE0=muE+JE*rE-JI*rI;
    II0=muI+JE*rE-JI*rI;
    mueffEvec(j)=IE0;
    mueffIvec(j)=II0;
    mueffEvec_ext(j+nextralowE)=IE0;
    mueffIvec_ext(j+nextralowI)=II0;
end

% Extend excitatory FI curve below attentional path
mueffE=mueffEvec(1);
for j=1:nextralowE
    rE=rEextendedvec(nextralowE-j+1);
    mueffE=fsolve(@(mueffE) fE(mueffE)-rE,mueffE,options);
    mueffEvec_ext(nextralowE-j+1)=mueffE;
end

% Extend inhibitory FI curve below attentional path
mueffI=mueffIvec(1);
for j=1:nextralowI
    rI=rIextendedvec(nextralowI-j+1);
    mueffI=fsolve(@(mueffI) fI(mueffI)-rI,mueffI,options);
    mueffIvec_ext(nextralowI-j+1)=mueffI;
end

% Extend excitatory FI curve above attentional path
histartE=length(rEvec)+nextralowE;
mueffE=mueffEvec(end);
for j=1:nextrahighE
    rE=rEextendedvec(histartE+j);
    mueffE=fsolve(@(mueffE) fE(mueffE)-rE,mueffE,options);
    mueffEvec_ext(histartE+j)=mueffE;
end

% Extend inhibitory FI curve above attentional path
histartI=length(rIvec)+nextralowI;
mueffI=mueffIvec(end);
for j=1:nextrahighI
    rI=rIextendedvec(histartI+j);
    mueffI=fsolve(@(mueffI) fI(mueffI)-rI,mueffI,options);
    mueffIvec_ext(histartI+j)=mueffI;
end

% plot excitatory FI curve with firing rate units of sp/s (multiply rates 
% by 100 to covert from units of sp/time constant, where time constant=1 
% is taken to be 10 ms)
figure; plot(mueffEvec_ext,100*rEextendedvec,'k')
hold on
% plot attentional path in red
plot(mueffEvec,100*rEvec,'r')
% mark beginning and end of attentional path as unattended and attended
% states
plot(mueffEvec(1),100*rEstart,'g.','MarkerSize',28)
plot(mueffEvec(end),100*rEend,'m.','MarkerSize',28)
hold off
xlim([0.3 0.8])
ylim([0 20])
xlabel('Input to exc. (units)')
ylabel('Exc. firing rate (sp/s)')
box off

% plot inhibitory FI curve with firing rate units of sp/s (multiply rates 
% by 100 to covert from units of sp/time constant, where time constant=1 
% is taken to be 10 ms)
figure; plot(mueffIvec_ext,100*rIextendedvec,'k')
hold on
% plot attentional path in red
plot(mueffIvec,100*rIvec,'r')
% mark beginning and end of attentional path as unattended and attended
% states
plot(mueffIvec(1),100*rIstart,'g.','MarkerSize',28)
plot(mueffIvec(end),100*rIend,'m.','MarkerSize',28)
hold off
xlim([0.3 0.8]);
ylim([0 25])
xlabel('Input to inh. (units)')
ylabel('Inh. firing rate (sp/s)')
box off

% plot the transformation from mu-eff to firing rates for E,I
figure; plot(mueffIvec,mueffEvec,'r','LineWidth',2)
xlabel(gca,'\mu_{eff,E}')
ylabel(gca,'\mu_{eff,I}')
hold on
plot(mueffIvec(1),mueffEvec(1),'g.','MarkerSize',28)
plot(mueffIvec(end),mueffEvec(end),'m.','MarkerSize',28)
xlim([0.4 0.65])
ylim([0.4 0.65])

% add axes on top and right
set(gca,'Box','off')
axesPosition=get(gca,'Position');
leftlim=get(gca,'ylim');
bottomlim=get(gca,'xlim');
A1=(rEvec(end)-rEvec(1))*(leftlim(1)-mueffEvec(1))/(mueffEvec(end)-mueffEvec(1))+rEvec(1);
A2=(leftlim(2)-mueffEvec(end))*(rEvec(end)-rEvec(1))/(mueffEvec(end)-mueffEvec(1))+rEvec(end);
B1=(rIvec(end)-rIvec(1))*(bottomlim(1)-mueffIvec(1))/(mueffIvec(end)-mueffIvec(1))+rIvec(1);
B2=(bottomlim(2)-mueffIvec(end))*(rIvec(end)-rIvec(1))/(mueffIvec(end)-mueffIvec(1))+rIvec(end);
hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YLim',[100*A1 100*A2],...            %#   ... and a different scale
                'YAxisLocation','right',...  %#   ... located on the righ
                'XLim',[100*B1 100*B2],...
                'XAxisLocation','top',...
                'Box','off');                %#   ... and no surrounding box
set(hNewAxes,'FontSize',14);
ylabel(hNewAxes,'Excitatory firing rate (Hz)');  %# Add a label to the right y axis
xlabel(hNewAxes,'Inhibitory firing rate (Hz)');
