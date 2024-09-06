%% File Name: Lab1.m
% Author: Abby Yoon   , Coby Garcia, Benji
% Created: September 15, 2023
% Description: Muscle contraction and activation lab
clear
close all
clc
%% Effect of Voluntary Change in Contractile Force (Figure 1)
abw=load('abby_weights.mat');
bn=1;
chn=3;
cn=2;
dt=1/abw.samplerate(1);
abwfrom=abw.datastart(chn,bn)+abw.com(cn,3);
abwto=abw.datastart(chn,bn)+abw.com(cn,3)+6*abw.samplerate(1)-1;
abwsnip=abw.data(abwfrom:abwto);
t=(0:length(abwsnip)-1)*dt;
chn=1;
abwfrom2=abw.datastart(chn,bn)+abw.com(cn,3);
abwto2=abw.datastart(chn,bn)+abw.com(cn,3)+6*abw.samplerate(1)-1;
abwsnip2=abw.data(abwfrom2:abwto2);
t=(0:length(abwsnip2)-1)*dt;
figure(1), plot(t,abwsnip)
hold on;
[yupper,ylower]=envelope(abwsnip,50,'rms');
plot(t,yupper);
plot(t,abwsnip2,'black')
xlabel('Time (s)');
ylabel('Normalized Amplitude (a.u.)');
legend('Raw EMG','ADI EMG','my RMS');
hold off;
%% Comparing RMS of Biceps Across All Group Members Figure 2
abw=load('abby_weights.mat');
anw=load('ana_weights.mat');
lyw=load('lydia_weights.mat');
bn=1;
chn=1;
cn=2;
dt=1/abw.samplerate(1);
abwfrom = abw.datastart(chn,bn);
abwto = abw.dataend(chn,bn);
abwsnip = abw.data(abwfrom:abwto);
t = [0:length(abwsnip)-1]*dt;
figure(2), plot(t,abwsnip)
hold on;
anwfrom = anw.datastart(chn,bn);
anwto = anw.dataend(chn,bn);
anwsnip = anw.data(anwfrom:anwto);
t = [0:length(anwsnip)-1]*dt;
plot(t,anwsnip)
hold on;
lywfrom = lyw.datastart(chn,bn);
lywto = lyw.dataend(chn,bn);
lywsnip = lyw.data(lywfrom:lywto);
t = [0:length(lywsnip)-1]*dt;
plot(t,lywsnip,'black')
xlabel('Time (s)');
ylabel('RMS Amplitude (mV)');
legend('Liv','Ana','Lydia');
hold off;
%% Compare Average RMS Biceps EMG Across Group Table 1 and Figure 3
abw=load('abby_weights.mat');
anw=load('ana_weights.mat');
lyw=load('lydia_weights.mat');
bn = 1;
chn = 1;
cn = 2;
dt = 1/abw.samplerate(1);
abwfrom=abw.datastart(chn,bn)+abw.com(cn,3);
abwto=abw.datastart(chn,bn)+abw.com(cn,3)+5*abw.samplerate(1)-1;
abwsnip=abw.data(abwfrom:abwto);
t=[0:length(abwsnip)-1]*dt;
abwrms1=rms(abwsnip);
cn=3;
abwfrom=abw.datastart(chn,bn)+abw.com(cn,3);
abwto=abw.datastart(chn,bn)+abw.com(cn,3)+5*abw.samplerate(1)-1;
abwsnip=abw.data(abwfrom:abwto);
t=[0:length(abwsnip)-1]*dt;
abwrms2=rms(abwsnip);
cn=4;
abwfrom=abw.datastart(chn,bn)+abw.com(cn,3);
abwto=abw.datastart(chn,bn)+abw.com(cn,3)+5*abw.samplerate(1)-1;
abwsnip=abw.data(abwfrom:liwto);
t=[0:length(abwsnip)-1]*dt;
abwrms3=rms(abwsnip);
cn=1;
dt = 1/anw.samplerate(1);
anwfrom=anw.datastart(chn,bn)+anw.com(cn,3);
anwto=anw.datastart(chn,bn)+anw.com(cn,3)+5*anw.samplerate(1)-1;
anwsnip=anw.data(anwfrom:anwto);
t=[0:length(anwsnip)-1]*dt;
anwrms1=rms(anwsnip);
cn=2;
anwfrom=anw.datastart(chn,bn)+anw.com(cn,3);
anwto=anw.datastart(chn,bn)+anw.com(cn,3)+5*anw.samplerate(1)-1;
anwsnip=anw.data(anwfrom:anwto);
t=[0:length(anwsnip)-1]*dt;
anwrms2=rms(anwsnip);
cn=3;
anwfrom=anw.datastart(chn,bn)+anw.com(cn,3);
anwto=anw.datastart(chn,bn)+anw.com(cn,3)+5*anw.samplerate(1)-1;
anwsnip=anw.data(anwfrom:anwto);
t=[0:length(anwsnip)-1]*dt;
anwrms3=rms(anwsnip);
cn=2;
dt = 1/lyw.samplerate(1);
lywfrom=lyw.datastart(chn,bn)+lyw.com(cn,3);
lywto=lyw.datastart(chn,bn)+lyw.com(cn,3)+5*lyw.samplerate(1)-1;
lywsnip=lyw.data(lywfrom:lywto);
t=[0:length(lywsnip)-1]*dt;
lywrms1=rms(lywsnip);
cn=3;
lywfrom=lyw.datastart(chn,bn)+lyw.com(cn,3);
lywto=lyw.datastart(chn,bn)+lyw.com(cn,3)+5*lyw.samplerate(1)-1;
lywsnip=lyw.data(lywfrom:lywto);
t=[0:length(lywsnip)-1]*dt;
lywrms2=rms(lywsnip);
cn=4;
lywfrom=lyw.datastart(chn,bn)+lyw.com(cn,3);
lywto=lyw.datastart(chn,bn)+lyw.com(cn,3)+5*lyw.samplerate(1)-1;
lywsnip=lyw.data(lywfrom:lywto);
t=[0:length(lywsnip)-1]*dt;
lywrms3=rms(lywsnip);
totlabrms=[abwrms1,abwrms2,abwrms3];
totanwrms=[anwrms1,anwrms2,anwrms3];
totlywrms=[lywrms1,lywrms2,lywrms3];
weight=[2.5,5,7.5];
figure(1),plot(weight,totabwrms,'*','LineStyle','--')
hold on;
plot(weight,totanwrms,'*','LineStyle','--')
plot(weight,totlywrms,'*','LineStyle','--')
ylim([0,5*10^-4])
xlabel('Weight (a.u.)');
ylabel('RMS Biceps Amplitude (mV)');
legend('Abby','Ana','Lydia');
hold off;
