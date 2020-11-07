function [ds ro4 bvec r tebal iteration]=dipoleinv(lap,f)
bvec1=f(:,1);
ds1=f(:,2);
bvec1(isnan(bvec1))=[];
ds1(isnan(ds1))=[];
bvec=bvec1';
ds=ds1';

for i=1:lap
    r(i)=mean(ds);
end

for i=1:lap-1
    tebal(i)=10;
end

ourvec=zeros(size(bvec));
diffs=zeros(size(bvec));
pm=[r tebal];
rawal=r;
tawal=tebal;
lr=length(r);
lt=length(tebal);
kr=10e-10;
iterasi=1;
maxiteration=500;
dfit=1;
iterasi=0;
iteration=0;

while iterasi<maxiteration
%forward
r1=forwarddip(r,tebal,bvec);
%inversi damped least square
beda=log(ds)-log(r1);
dd=beda;
misfit1=beda*beda';
    if misfit1<kr
        loglog(bvec,ds,'ro',bvec,r1,'b');
        axis([1 1000 1 1000])
        xlabel('AB/2(m)');
        ylabel('Apparent Resistivity (Ohm-m)');
        break
    end
%jacobian
    [A]=jacobian(r,tebal,lr,lt,ds,r1,bvec);
    [U S V]=svd(A,0);
    ss=length(S);
    say=1;
    k=0;

    while say<ss
        diagS=diag(S);
        beta=diagS(say)*(dfit^(1/say));
        if beta<10e-5
            beta=0.001*say;
        end
        for i4=1:ss
                SS(i4,i4)=S(i4,i4)/(S(i4,i4)^2+beta);
        end
        dmg=V*SS*U'*dd';
        mg=exp(log(pm)+dmg');
        r=mg(1:lr);
        tebal=mg(1+lr:lr+lt);
        ro4=forwarddip(r,tebal,bvec);
        beda2=log(ds)-log(ro4);
        misfit2=beda2*beda2';
        if misfit2>misfit1 
            say=say+1
            k=k+1;
            if k==ss-1
                iterasi=maxiteration; 
                say=ss+1;
                break
            end
            else
                say=ss+1;
                pm=mg;
                dfit=(misfit1-misfit2)/misfit1;
                iterasi=iterasi+1;
                a=iterasi;
                if dfit<kr
                    iterasi=maxiteration;
                    say=say+1;
                    b=b+1;
                end
        end
    end
iteration=iteration+1;
end
end
