function [ nrm,grad1 ] = covgivenpars_of( g,eqind,N )

nn=length(eqind);

grad1=zeros(N,1);
CC=zeros(nn,1);
for i=1:nn
    x=eqind(i,1);
    y=eqind(i,2);
    cu=eqind(i,3);
    ca=eqind(i,4);
    gi=g(x); gj=g(y);
    base=gi*gj*cu-ca;
    CC(i)=base;
    grad1(x)=grad1(x)+gj*cu*base;
    grad1(y)=grad1(y)+gi*cu*base;
end

nrm=norm(CC);
grad1=grad1/nrm;

end