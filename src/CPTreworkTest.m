clc
clear
close all

Fs = 10e3;
Ts = 1/Fs;
T_max = 0.3;
t= 0:Ts:T_max-Ts;                   % 10 seconds @ 1kHz sample rate
f0 = 5;
f1 = 150;
a = linspace(2,0.5,length(t));
% y = a.*chirp(t,f0,T_max,f1) + 0.0*randn(1,length(t));

a1 = 1.5; a2 = 1;
y = a1*cos(2*pi*10*t) + a2*cos(2*pi*60*t) + 0.0000*randn(1,length(t));
Hy = a1*sin(2*pi*10*t) + a2*sin(2*pi*60*t);
f_max = 60;
figure
plot(y,Hy,'-+')
drawnow

figure
plot(y)

InitialPoints = 4;
PointsStep = 1;
UpperLimit = 10;
PlotMode = 1;

[IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0, m_star, M] = cpt(y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode,f_max);
% N_imfs = 3;
% [C r_approx IF_interp phi_interp phi_unwrapped m_star, M] = CPT_EMD_rework(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);

% figure
% plot(t,C,t,y)
% 
% figure
% plot(t,sum(C,1))

figure
ax(1) = subplot(211);
plot(t,x,t,y)
legend('x','y')

ax(2) = subplot(212);
plot(t,M,'.')

linkaxes(ax,'x')
