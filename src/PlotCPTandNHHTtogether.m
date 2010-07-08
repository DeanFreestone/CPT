clear
clc
close all

% FS = 18;                % fontsize for figures
FS2 = 10;                   % for the legend
FS3 = 25;

Fs = 10e3;
Ts = 1/Fs;
T_max = 0.3;
t= 0:Ts:T_max-Ts;                   % 10 seconds @ 1kHz sample rate

f0 = 10;
f1=100;                    % Start at 1 Hz, go up to 150 Hz

f_middle = f0+(f1-f0)/2;
f = linspace(f0,f1,length(t));
% true_IF_middle = f(Sample2Cut+1:end-Sample2Cut);

true_IP = zeros(1,length(t));
for n=2:length(t)
    true_IP(n) = true_IP(n-1) + f(n)*2*pi*Ts;
end

true_IP = mod(true_IP,2*pi);

A1 = 4;
a = linspace(A1,.5,length(t));
y = a.*chirp(t,f0,T_max,f1);
% y_middle = y(Sample2Cut+1:end-Sample2Cut);

sigma = 0.00:0.07:0.6;

PointMultiplier = 5;

MinPoints = round(0.9*Fs/f1);
MaxPoints = round(PointMultiplier*MinPoints);
PointsStep = 2;
PlotMode = 0;

PhaseFig  = figure('units','normalized','position',[0 0 0.4 0.6]);
FS = 16;
LW = 1.6;
MS = 5;
MS2 = 3;

EdgeEffectTime = 0.02;
EdgeEffectSample = round(EdgeEffectTime*Fs);
t_start_sample = 1;
NIterations = 100;

% snr_db_start = zeros(NIterations,length(sigma));
% snr_db_mid = zeros(NIterations,length(sigma));
% snr_db_end= zeros(NIterations,length(sigma));
% RMSE_start = zeros(NIterations,length(sigma));
% RMSE_mid = zeros(NIterations,length(sigma));
% RMSE_end = zeros(NIterations,length(sigma));

s_1 = round(Fs*T_max/3);
s_2 = round(2*Fs*T_max/3);

plotmode = 0;
tic
n=3;
noise = sigma(n)*randn(1,length(t));
y_n = y + noise;

rms_y_start = norm(y(1:s_1))/sqrt(length(y(1:s_1)));
rms_n_start = norm(noise(1:s_1))/sqrt(length(noise(1:s_1)));

rms_y_mid = norm(y(s_1:s_2))/sqrt(length(y(s_1:s_2)));
rms_n_mid = norm(noise(s_1:s_2))/sqrt(length(noise(s_1:s_2)));

rms_y_end = norm(y(s_2:end))/sqrt(length(y(s_2:end)));
rms_n_end = norm(noise(s_2:end))/sqrt(length(noise(s_2:end)));

snr_start = (rms_y_start/rms_n_start)^2;
snr_db_start = 10*log10(snr_start);

snr_mid = (rms_y_mid/rms_n_mid)^2;
snr_db_mid = 10*log10(snr_mid);

snr_end = (rms_y_end/rms_n_end)^2;
snr_db_end = 10*log10(snr_end);

[IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0, m_star, M] = cpt_rework(y_n, MinPoints, MaxPoints, PointsStep, Ts, PlotMode);

phi_est_cpt = mod(phi_interp,2*pi);

phi_diff = min( (repmat(true_IP,3,1) - [phi_est_cpt-2*pi ; phi_est_cpt ; phi_est_cpt + 2*pi]).^2 ,[],1); 
phi_diff_start = phi_diff(EdgeEffectSample:s_1);
phi_diff_mid = phi_diff(s_1:s_2);
phi_diff_end = phi_diff(s_2:end-EdgeEffectSample);

RMSE_start_cpt = sqrt(mean(phi_diff_start));
RMSE_mid_cpt = sqrt(mean(phi_diff_mid));
RMSE_end_cpt = sqrt(mean(phi_diff_end));

% ~~~~~
ax(1) = subplot(3,1,1,'parent',PhaseFig);
plot(t,y_n,'k','parent',ax(1))
grid(ax(1),'on')
Pos = get(ax(1),'position');
% set(ax(1),'position',[Pos(1) Pos(2) Pos(3) Pos(4)]);
title(ax(1),['Mean SNRs = ' num2str(snr_db_start) ', ' num2str(snr_db_mid) ', ' num2str(snr_db_end) ' dB'],'fontsize',FS)%'interpreter','latex')
ylabel(ax(1),'Amplitude','fontsize',FS,'position',[-0.031 -0.5 1])%,'interpreter','latex')
set(ax(1),'xticklabel',[],'fontsize',FS,'ytick',[-2 0 2 4])
xlim(ax(1),[t(t_start_sample) t(end)])
ylim(ax(1),[-A1 A1])

ax(2) = subplot(3,1,2,'parent',PhaseFig);
plot(t,true_IP,'k.','linewidth',LW,'markersize',MS,'parent',ax(2))
hold(ax(2),'on')
plot(t,phi_est_cpt,'.r','markersize',MS,'parent',ax(2))
hold(ax(2),'off')
Pos = get(ax(2),'position');
% set(ax(2),'position',[Pos(1) Pos(2)+0.01 Pos(3) Pos(4)+0.02]);
grid(ax(2),'on')
set(ax(2),'xticklabel',[])
set(ax(2),'fontsize',FS,'ytick',[0 pi 2*pi],'yticklabel',{'0', '3.1','6.2'})
ylabel(ax(2),'IP CPT (rads)','fontsize',FS,'position',[-0.031 3 1])%,'interpreter','latex')
xlim(ax(2),[t(t_start_sample) t(end)])
ylim(ax(2),[0,2*pi+0.4])

% ~~~~~~
% ~~~~~~

imf = emd_orig(y_n);
for nnn=1:size(imf,2)
            
    imf_HHT = imf{nnn};                  % need to normalise these guys    
    NormIMF = NormalizeIMFs(imf_HHT);
    phi_NHHT = atan2(imag(hilbert(NormIMF)),NormIMF);
    phi_est_NHHT(nnn,:) = mod(phi_NHHT,2*pi);

    phi_diff = min( (repmat(true_IP,3,1) - [phi_est_NHHT(nnn,:)-2*pi ; phi_est_NHHT(nnn,:) ; phi_est_NHHT(nnn,:) + 2*pi]).^2 ,[],1); 
    phi_diff_start = phi_diff(EdgeEffectSample:s_1);
    phi_diff_mid = phi_diff(s_1:s_2);
    phi_diff_end = phi_diff(s_2:end-EdgeEffectSample);

    RMSE_start_all(nnn) = sqrt(mean(phi_diff_start));
    RMSE_mid_all(nnn) = sqrt(mean(phi_diff_mid));
    RMSE_end_all(nnn) = sqrt(mean(phi_diff_end));

end

[RMSE_mid_NHHT I] = min(RMSE_mid_all);
RMSE_start_NHHT = RMSE_start_all(I);
RMSE_end_NHHT = RMSE_end_all(I);
phi_est_NHHT_best = phi_est_NHHT(I,:);

ax(3) = subplot(3,1,3,'parent',PhaseFig);
plot(t,true_IP,'k.','linewidth',LW,'markersize',MS,'parent',ax(3))
hold(ax(3),'on')
plot(t,phi_est_NHHT_best,'.r','markersize',MS,'parent',ax(3))
hold(ax(3),'off')
Pos = get(ax(3),'position');
% set(ax(2),'position',[Pos(1) Pos(2)+0.01 Pos(3) Pos(4)+0.02]);
grid(ax(3),'on')
set(ax(3),'fontsize',FS,'ytick',[0 pi 2*pi],'yticklabel',{'0', '3.1','6.2'})
ylabel(ax(3),'IP NHHT (rads)','fontsize',FS,'position',[-0.031 3 1])%,'interpreter','latex')
xlabel(ax(3),'Times (s)','fontsize',FS)%,'interpreter','latex')
xlim(ax(3),[t(t_start_sample) t(end)])
ylim(ax(3),[0,2*pi+0.4])