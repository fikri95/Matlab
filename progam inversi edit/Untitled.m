close all
clear all 
clc

lap=input('jumlah lapisan :');
[ds ro4 bvec r tebal iteration]=tschrinv(lap);
subplot(1,2,1),
loglog(bvec,ds,'ro',bvec,ro4,'b');
axis([1 1000 1 1000])
xlabel('AB/2(m)');
ylabel('Apparent Resistivity (Ohm-m)');


observed=ds;
calculated=ro4;
legend('obs','cal');
format bank;
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