cvalvec=zeros(ndraws,1);
cprevec=zeros(ndraws,1);

for ix=1:ndraws
    eqrem=squeeze(eqremvec(ix,:));
    g=squeeze(gvecall(ix,:));
    xind=eqrem(1);
    yind=eqrem(2);
    cu=eqrem(3);
    ca=eqrem(4);
    cpre=g(xind)*g(yind)*cu;
    cvalvec(ix)=ca;
    cprevec(ix)=cpre;     
end

cvalbest=cvalvec(:,1);
cprebest=cprevec(:,1);

r2best=corrcoef(cprebest,cvalbest);
%{
CPB=cprebest; CVB=cvalbest;



P=(CPB-CVB).^2;
CPV=cprevec;
CVV=cvalvec;

figure; hold on
for i=1:length(CPB)
    plot(CPV(i,2),P(i),'.');
    plot(CPV(i,3),P(i),'.');
    plot([CPV(i,2);CPV(i,3)],[P(i);P(i)],'r')
end

box=zeros(N,2*(N-1));
boxi=zeros(N,1)+1;
for i=1:length(CPB)
    n=CPV(i,2);
    box(n,boxi(n))=P(i);
    boxi(n)=boxi(n)+1;
    n=CPV(i,3);
    box(n,boxi(n))=P(i);
    boxi(n)=boxi(n)+1;
end
boxm=mean(box,2);
figure; plot(1:15,boxm,'.')


figure; hold on
set(gca,'FontSize',18)
%scatter(cpremeanU,cvalmeanU,'bo');
plot(cprebest,cvalbest,'b.','MarkerSize',24);
%scatter(cpremeanA,cvalmeanA,'go');
%x=-50:0.1:350;
x=-100:0.1:500;
plot(x,x,'k','LineWidth',2)
hold off
xlabel('Predicted Covariance')
ylabel('Actual Covariance')
%legend('Unattended','Attended')
%xlim([-100 500])
%}