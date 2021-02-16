N=length(covmatatt);

scma=sqrtm(covmatatt);
cap=[];
for i=1:N
    for j=i+1:N
        cap=[cap scma(i,j)];
    end
end       
qa=cap(randperm(length(cap)));
qadd=diag(scma);
qad=qadd(randperm(length(qadd)));
qqa=zeros(N);
ct=1;
for i=1:N
    for j=1:N
        if i<j
            qqa(i,j)=qa(ct);
            qqa(j,i)=qa(ct);
            ct=ct+1;
        elseif i==j
            qqa(i,j)=qad(i);
        end
    end
end
covmatatt=real(qqa*qqa);
