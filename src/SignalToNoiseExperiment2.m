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

f_max = f1;

f_middle = f0+(f1-f0)/2;
f = linspace(f0,f1,length(t));
% true_IF_middle = f(Sample2Cut+1:end-Sample2Cut);

true_IP = zeros(1,length(t));
for n=2:length(t)
    true_IP(n) = true_IP(n-1) + f(n)*2*pi*Ts;
end

true_IP = mod(true_IP,2*pi);

A1 = 4;
A2 = 1;
Envelope = linspace(A1,A2,length(t));
% Envelope2 = [Envelope(1:round(length(Envelope)/2)) fliplr(Envelope(1:round(length(Envelope)/2)))];
Envelope2 = [linspace(A1,A2,round(length(Envelope)/2)) linspace(A2,A1,round(length(Envelope)/2))];

y = Envelope2.*chirp(t,f0,T_max,f1);
figure
plot(t,y,t,Envelope)
% y_middle = y(Sample2Cut+1:end-Sample2Cut);

sigma = 0.0:0.07:0.6;

PointMultiplier = 5;
MinPoints = round(0.9*Fs/f1);

MaxPoints = round(PointMultiplier*MinPoints);
PointsStep = 5;
PlotMode = 0;

MaxPoints = floor((1/f_middle)*Fs/2)+100;
MinPoints = floor((1/f_middle)*Fs/2);
PointsStep = 20;

FS = 16;
LW = 1.6;
MS = 5;

EdgeEffectTime = 0.02;
EdgeEffectSample = round(EdgeEffectTime*Fs);
t_start_sample = 1;
NIterations = 10;

f_max = f(end-EdgeEffectSample);

snr_db_start = zeros(NIterations,length(sigma));
snr_db_mid = zeros(NIterations,length(sigma));
snr_db_end= zeros(NIterations,length(sigma));
RMSE_start = zeros(NIterations,length(sigma));
RMSE_mid = zeros(NIterations,length(sigma));
RMSE_end = zeros(NIterations,length(sigma));

s_1 = round(Fs*T_max/3);
s_2 = round(2*Fs*T_max/3);

plotmode = 1;
tic
PhaseFig  = figure('units','normalized','position',[0 0 1 1]);

for n=1:length(sigma)
    disp(['iteration ' num2str(n) ' of ' num2str(length(sigma)) ' noise levels']) 
    
    for nn=1:NIterations
        
        noise = sigma(n)*randn(1,length(t));
        y_n = y + noise;

        rms_y_start = norm(y(1:s_1))/sqrt(length(y(1:s_1)));
        rms_n_start = norm(noise(1:s_1))/sqrt(length(noise(1:s_1)));

        rms_y_mid = norm(y(s_1:s_2))/sqrt(length(y(s_1:s_2)));
        rms_n_mid = norm(noise(s_1:s_2))/sqrt(length(noise(s_1:s_2)));

        rms_y_end = norm(y(s_2:end))/sqrt(length(y(s_2:end)));
        rms_n_end = norm(noise(s_2:end))/sqrt(length(noise(s_2:end)));

        snr_start = (rms_y_start/rms_n_start)^2;
        snr_db_start(nn,n) = 10*log10(snr_start);

        snr_mid = (rms_y_mid/rms_n_mid)^2;
        snr_db_mid(nn,n) = 10*log10(snr_mid);

        snr_end = (rms_y_end/rms_n_end)^2;
        snr_db_end(nn,n) = 10*log10(snr_end);

%         [IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0, m_star, M] ...
%             = cpt(y_n, MinPoints, MaxPoints, PointsStep, Ts, PlotMode, f_max);
%         phi_est = mod(phi_interp,2*pi);

        DesiredArcLength = 3*pi/4;
        theta_diff_thresh = pi;
        [phi r] = CPTfunction(y_n,Ts,DesiredArcLength,theta_diff_thresh,f_max);
        phi_est = mod(phi, 2*pi);
        
        phi_diff = min( (repmat(true_IP,3,1) - [phi_est-2*pi ; phi_est ; phi_est + 2*pi]).^2 ,[],1); 
        phi_diff_start = phi_diff(EdgeEffectSample:s_1);
        phi_diff_mid = phi_diff(s_1:s_2);
        phi_diff_end = phi_diff(s_2:end-EdgeEffectSample);

        RMSE_start(nn,n) = sqrt(mean(phi_diff_start));
        RMSE_mid(nn,n) = sqrt(mean(phi_diff_mid));
        RMSE_end(nn,n) = sqrt(mean(phi_diff_end));

        if plotmode == 1
            ax(n+floor((n-1)/3)*3) = subplot(6,3,n+floor((n-1)/3)*3,'parent',PhaseFig);
            plot(t,y_n,'k','parent',ax(n+floor((n-1)/3)*3))
            grid(ax(n+floor((n-1)/3)*3),'on')
            Pos = get(ax(n+floor((n-1)/3)*3),'position');
            set(ax(n+floor((n-1)/3)*3),'position',[Pos(1) Pos(2) Pos(3) Pos(4)-0.01]);
            title(ax(n+floor((n-1)/3)*3),['SNR = ' num2str(snr_db_start(nn,n)) ', ' num2str(snr_db_mid(nn,n)) ', ' num2str(snr_db_end(nn,n)) ' dB'],'fontsize',FS,'fontname','Times New Roman')%'interpreter','latex')
            if mod(n,3) ~= 1
                set(ax(n+floor((n-1)/3)*3),'yticklabel',[])
            else
                ylabel(ax(n+floor((n-1)/3)*3),'Amplitude','fontsize',FS,'fontname','Times New Roman','position',[-0.031 -0.5 1])%,'interpreter','latex')
            end
            set(ax(n+floor((n-1)/3)*3),'xticklabel',[],'fontsize',FS,'fontname','Times New Roman','ytick',[-2 0 2 4])
            xlim(ax(n+floor((n-1)/3)*3),[t(t_start_sample) t(end)])
            ylim(ax(n+floor((n-1)/3)*3),[-A1 A1])

            ax(6,3,n+(floor((n-1)/3)+1)) = subplot(6,3,n+(floor((n-1)/3)+1)*3,'parent',PhaseFig);
            plot(t,true_IP,'k','linewidth',LW,'parent',ax(6,3,n+(floor((n-1)/3)+1)))
            hold(ax(6,3,n+(floor((n-1)/3)+1)),'on')
            plot(t,phi_est,'.r','markersize',MS,'parent',ax(6,3,n+(floor((n-1)/3)+1)))
            hold(ax(6,3,n+(floor((n-1)/3)+1)),'off')
            Pos = get(ax(6,3,n+(floor((n-1)/3)+1)),'position');
            set(ax(6,3,n+(floor((n-1)/3)+1)),'position',[Pos(1) Pos(2)+0.01 Pos(3) Pos(4)+0.02]);
            grid(ax(6,3,n+(floor((n-1)/3)+1)),'on')
            if n<7
                set(ax(6,3,n+(floor((n-1)/3)+1)),'xticklabel',[])
            else
                xlabel(ax(6,3,n+(floor((n-1)/3)+1)),'Times (s)','fontsize',FS,'fontname','Times New Roman')%,'interpreter','latex')
            end
            set(ax(6,3,n+(floor((n-1)/3)+1)),'fontsize',FS,'fontname','Times New Roman','ytick',[0 pi 2*pi],'yticklabel',{'0', '3.14','6.28'})
            if mod(n,3) ~= 1
                set(ax(6,3,n+(floor((n-1)/3)+1)),'yticklabel',[])
            else
                ylabel(ax(6,3,n+(floor((n-1)/3)+1)),'Radians','fontsize',FS,'fontname','Times New Roman','position',[-0.031 3 1])%,'interpreter','latex')
            end
            xlim(ax(6,3,n+(floor((n-1)/3)+1)),[t(t_start_sample) t(end)])
            ylim(ax(6,3,n+(floor((n-1)/3)+1)),[0,2*pi+0.4])
            drawnow
        end
    end
end
toc

y_max = 2;
FS = 16;
roundlabelstart = round(10*mean(snr_db_start,1))/10;
roundlabelmid = round(10*mean(snr_db_mid,1))/10;
roundlabelend = round(10*mean(snr_db_end,1))/10;
for n=1:length(sigma)
    xticklabelstart{n} = num2str(roundlabelstart(n));
    xticklabelmid{n} = num2str(roundlabelmid(n));
    xticklabelend{n} = num2str(roundlabelend(n));
end
xticklabelstart = fliplr(xticklabelstart);
xticklabelmid = fliplr(xticklabelmid);
xticklabelend = fliplr(xticklabelend);
y_max = 50;
subplotextension = 0.055;
boxplotfig = figure('units','normalized','position',[0 0 1 1]);
ax(1) = subplot(1,3,1,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_start*100/pi,roundlabelstart,'parent',ax(1))
ylim(ax(1),[0 y_max])
pos = get(ax(1),'position');
set(ax(1),'xtick',1:9,'xticklabel',xticklabelstart,'fontsize',FS,'position',[pos(1)-subplotextension/2 pos(2) pos(3)+subplotextension pos(4)])
ylabel(ax(1),'RMSE CPT','fontsize',FS)

ax(2) = subplot(1,3,2,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_mid*100/pi,roundlabelmid,'parent',ax(2))
ylim(ax(2),[0 y_max])
pos = get(ax(2),'position');
set(ax(2),'xtick',1:9,'xticklabel',xticklabelmid,'yticklabel',{''},'fontsize',FS,'position',[pos(1)-subplotextension/2 pos(2) pos(3)+subplotextension pos(4)])

ax(3) = subplot(1,3,3,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_end*100/pi,roundlabelend,'parent',ax(3))
ylim(ax(3),[0 y_max])
pos = get(ax(3),'position');
set(ax(3),'xtick',1:9,'xticklabel',xticklabelend,'yticklabel',{''},'fontsize',FS,'position',[pos(1)-subplotextension/2 pos(2) pos(3)+subplotextension pos(4)])

