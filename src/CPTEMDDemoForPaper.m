clear 
clc
close all

Fs = 5e3;                           % samples/second
Ts = 1/Fs;                          % sample period (radians)
Duration = 3;                       % seconds
NSamples = Duration*Fs;
t = linspace(0,Duration,NSamples);

f1 = 7;                    Theta1 = 2*pi*f1*t;                % frequency Hz
f2 = 13;                  Theta2 = 2*pi*f2*t;

alpha = 1;
A1 = 1/(f1^alpha);              % power falls off at 1/(f^2) and amplitude falls away at 1/f
A2 = 1/(f2^alpha);   

x1 = A1*cos(Theta1);
x2 = A2*cos(Theta2);

y = x1 + x2;

InitialPoints = 5;
UpperLimit = 40;
PointsStep = 1;
N_imfs = 1;
PlotMode = 0;

[C r_approx IF_interp phi_interp phi_unwrapped m_star, M] = ...
    CPT_EMD_rework(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);