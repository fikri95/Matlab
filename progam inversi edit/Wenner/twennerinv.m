function [ds ro4 bvec r tebal iteration]=twennerinv(lap,f)
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
%s=2:10:1000;
%bvec=s./2;


ourvec=zeros(size(bvec));
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
%data sintetik
%ds=datasintetik(s,bvec,ruji,tuji);
   
%yd=y+rand(1,numel(y)).*y/100     
%ds1=forwardw(ruji,tuji,bvec);
%pest=ds1*10;
    %ds=ds1+rand(1,numel(ds1)).*pest/100;
while iterasi<maxiteration
%forward
r1=forwardw(r,tebal,bvec);
%inversi damped least square
beda=[log(ds)-log(r1)];
dd=beda;
misfit1=beda*beda';
    if misfit1<kr
        ro4=r1;
        break
    end
%jacobian
    [A]=jacobianw(r,tebal,lr,lt,ds,r1,bvec);
    [U S V]=svd(A,0);
    ss=length(S);
    el=1;
    k=0;

    while el<ss
        diagS=diag(S);
        epsilon=diagS(el)*(dfit^(1/el)); %modifikasi epsilon
        if epsilon<10e-5
            epsilon=0.001*el;
        end
        for i4=1:ss
                SS(i4,i4)=S(i4,i4)/(S(i4,i4)^2+epsilon);
        end
        dmg=V*SS*U'*dd';
        modelhasil=exp(log(pm)+dmg');
        r=modelhasil(1:lr);
        tebal=modelhasil(1+lr:lr+lt);
        ro4=forwardw(r,tebal,bvec);
        beda2=[log(ds)-log(ro4)];
        misfit2=beda2*beda2';
        if misfit2>misfit1 
            el=el+1;
            k=k+1;
            if k==ss-1
                iterasi=maxiteration; 
                el=ss+1;
              
                break
            end
            else
                el=ss+1;
                pm=modelhasil;
                dfit=(misfit1-misfit2)/misfit1;
                iterasi=iterasi+1;
                a=iterasi;
                if dfit<kr
                    iterasi=maxiteration;
                    el=el+1;
                  
                end
        end
    end

iteration=iteration+1;
end
end