close all
clear all 
clc
%input parameter
lap=input('jumlah lapisan :');
for i=1:lap
    r(i)=input(['resistivitas lapisan ke-',num2str(i),':']);
end

for i=1:lap-1
    tebal(i)=input(['tebal lapisan ke-',num2str(i),':']);
end
%s=2:10:1000;
%bvec=s./2;
bvec1=xlsread('Book1','A:A');
ds1=xlsread('Book1','C:C');
bvec=bvec1';
ds=ds1';
%ruji=[250 80 100 200];
%tuji=[10 20 20];
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
%data sintetik
%ds=datasintetik(s,bvec,ruji,tuji);
%ds=forward(ruji,tuji,bvec);
%ds1=forward(ruji,tuji,bvec);
%pest=ds1*10;
    %ds=ds1+rand(1,numel(ds1)).*pest/100;
while iterasi<maxiteration
%forward
r1=forward(r,tebal,bvec);
%inversi damped least square
beda=[log(ds)-log(r1)];
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
        ro4=forward(r,tebal,bvec);
        beda2=[log(ds)-log(ro4)];
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
subplot(1,2,1),
loglog(bvec,ds,'ro',bvec,ro4,'b');
axis([1 1000 1 1000])
xlabel('AB/2(m)');
ylabel('Apparent Resistivity (Ohm-m)');
pause(0.001)
iteration=iteration+1;
end
observed=ds;
calculated=ro4;
legend('obs','cal');
format bank;
pm
rr=[0,r];
tt=[0,cumsum(tebal),max(tebal)*10];
subplot(1,2,2),
stairs(rr,tt,'r-');
%rrr=[0,ruji];
%ttt=[0,cumsum(tuji),max(tebal)*10];
hold on;
subplot(1,2,2),
%stairs(rrr,ttt,'b--');
set(gca,'Ydir','reverse');
set(gca,'Xscale','log');
ylim([0 150]);
xlim([0 1000]);
xlabel('Resistivity (Ohm-m)');
ylabel('Depth (m)');
rms=norm(ro4-ds)/sqrt(length(ds));
title (['\bf \fontsize{14}\fontname{Times}iterasi = ',...
    num2str(iteration),' ; rms = ', num2str(rms)]);