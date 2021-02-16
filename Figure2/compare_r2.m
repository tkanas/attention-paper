cvalvec=zeros(ncv,1);
cprevec=zeros(ncv,1);
for ix=1:ncv
    xind=eqind(ix,1);
    yind=eqind(ix,2);
    cu=eqind(ix,3);
    ca=eqind(ix,4);
    cpre=g(xind)*g(yind)*cu;
    cvalvec(ix)=ca;
    cprevec(ix)=cpre;        
end
r2best=corrcoef(cprevec,cvalvec);