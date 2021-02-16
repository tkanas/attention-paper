N=length(covmatatt);
ncv=N*(N-1)/2;

xvec=zeros(N,1);
fvec=0;
gvecall=zeros(ndraws,N);
fvecall=zeros(ndraws,1);

%eqindvec=zeros(ndraws,nknown/2,4);
eqremvec=zeros(ndraws,4);

eqall=[];
for i=1:N
    for j=i+1:N
        eqall=[eqall;[i,j,covmatunatt(i,j),covmatatt(i,j)]];
    end
end

rp=randperm(ncv,ndraws);
for i=1:ndraws
    %if rem(i,100)==0
    %    fprintf('system %d out of %d \n',i,ndraws)
    %end
    eqind=eqall;
    eqrem=squeeze(eqall(rp(i),:));
    eqind(rp(i),:)=[];
    
    %eqindvec(i,1:size(eqind,1),:)=eqind;
    eqremvec(i,:)=eqrem;
    
    input0=rand(N,1); % gUN,gAT,Srel

    covanonf_of=@(inputs)covgivenpars_of(inputs,eqind,N);
    opts=optimset('GradObj','on','Display','off');

    [g,f] = fminunc(covanonf_of,input0,opts);
        
    gvecall(i,:)=g;
    fvecall(i)=f;
    
end