
% plot eeg for paper

clc
clear
close all

% load CSGridDepth200608CS33_19900_Chs_1_48_200s
load HGSP7_98s_seizure
Fs = 4069.010498046875;     % Hz


% pick channel
x = detrend(Data(1:end-round(Fs*2),end) - Data(1:end-round(2*Fs),end-1));
% x = detrend(Data(:,end-1));

t = 0:1/Fs:(length(x)-1)/Fs;        % seconds

clear data

% use a median filter to give the data a first clean
x = medfilt1(x,20);

% low pass filter file at 40 Hz
Fc = [2.5 95];                                % Hz
Wc = Fc/(Fs/2);                     % normalised digital frequency
[b a] = butter(3,Wc);
y = filtfilt(b,a,x);

figure
ax(1) = subplot(211);
plot(t,x)
ax(2) = subplot(212);
plot(t,y)
linkaxes(ax,'x')

% parameters for the cpt
N_imfs = 4;
InitialPoints = 6;                 %~22 samples is half the period of 90 Hz
UpperLimit = 40;
PointsStep = 2;
Ts = 1/Fs;
PlotMode = 0;

figure
tic
[C r_approx IF_interp phi_interp phi_unwrapped m_star, M] =...
    CPT_EMD_rework(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);
toc

figure
ax1(1) = subplot(N_imfs+1,1,1);
plot(t,detrend(y))
for n=1:N_imfs
    ax1(n+1) = subplot(N_imfs+1,1,n+1);
    plot(t,C(n,:))
end
linkaxes(ax1,'x')

diff_phi = diff(phi_unwrapped');
IF = [zeros(1,N_imfs); diff_phi]/(2*pi*Ts);
figure
plot(t,IF(:,1:N_imfs),'.')



%%
FS = 20;
ImageFig = figure('units','normalized','position',[0 0 1 1]);
ax1(1) = subplot(5,1,1,'parent',ImageFig);
plot(t,y,'k')
xlim(ax1(1),[t(1) t(end)])
ylim(ax1(1),[min(y) max(y)])
set(ax1(1),'xticklabel',[],'fontsize',FS,'ytick',[-1e-3 0 1e-3 2e-3],'yticklabel',{'-1','0','1','2'})
ylabel(ax1(1),'iEEG (mV)','fontsize',FS,'position',[-7 0 1])

TimeSampleStart = 6;
TimeSampleEnd = 6;
for n=1:N_imfs
    [t_out Freq_Vector IF_image IA_IF_image] = STIFH(IF_interp(n,:),r_approx(n,:),Fs,200,0,100,0.15,75);

    ax1(n+1) = subplot(N_imfs+1,1,n+1,'parent',ImageFig);
    imagesc(t_out(TimeSampleStart:end-TimeSampleEnd),...
    Freq_Vector(1:end-1),IF_image(1:end-1,TimeSampleStart:end-TimeSampleEnd),'parent',ax1(n+1))
    set(ax1(n+1),'ydir','normal')
    if n<N_imfs
        set(ax1(n+1),'xticklabel',[])
    else
        xlabel(ax1(n+1),'Time (s)','fontsize',FS)
        ylabel(ax1(n+1),'Frequency (Hz)','fontsize',FS,'position',[-7 15.23 9.16])
    end
    
    Pos = get(ax1(n+1),'position');
    set(ax1(n+1),'position',[Pos(1) Pos(2) Pos(3) Pos(4)+0.04],'fontsize',FS)
    
    colormap(ax1(n+1),'hot')
    colorbar('position',[Pos(1)+Pos(3)+0.01 Pos(2) 0.02 Pos(4)+0.02],'fontsize',FS)
end

ylim(ax1(2),[0 100]), ylim(ax1(3),[0 50]), ylim(ax1(4),[0 25]),ylim(ax1(5),[0 12.5])
set(ax1(2),'ytick',[20 40 60],'yticklabel',{'20','40','60'})
set(ax1(3),'ytick',[10 20 30])
% set(ax1(4),'ytick',[5 10 15])
% set(ax1(5),'ytick',[2.5 5 7.5])
linkaxes(ax1,'x')

y=y*1e4;
C = C*1e4;

% ~~~~~~~~~~~
%%
TimeSeriesFig = figure('units','normalized','position',[0 0 1 1]);
NRows = 5;
NCols = 4;

m=1;
ax2(m) = subplot(NRows,1,1,'parent',TimeSeriesFig);
plot(t,y,'k')
xlim(ax2(m),[t(1) t(end)])
ylim(ax2(m),[min(y) max(y)])
set(ax2(m),'fontsize',FS,'ytick',[-5 0 5],'yticklabel',{'-5','0','5'})
WindowSize = 0.5;

% plot the first feature
t_start = 1.2;
t_end = t_start + WindowSize;
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'g','edgecolor','none','parent',ax2(1),'facealpha',0.4)

s_start = round(t_start*Fs);
s_end = round(t_end*Fs);

m=m+1;
ax2(m) = subplot(NRows,NCols,5,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),y(s_start:s_end),'k','parent',ax2(m))
hold(ax2(m),'on')
plot(t(s_start:s_end),sum(C(:,s_start:s_end),1),'r','parent',ax2(m))
hold(ax2(m),'off')
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'ytick',[-5 0 5],'yticklabel',{'-5','0','5'})
text(1.2, 0, '$y(n)$','interpreter','latex','fontsize',FS)
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'g','edgecolor','none','parent',ax2(m),'facealpha',0.4)

m=m+1;
ax2(m) = subplot(NRows,NCols,9,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(1,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'ytick',[-5 0 5],'yticklabel',{'-5','0','5'})
text(1.2, 0, '$x_1(n)$','interpreter','latex','fontsize',FS)

m=m+1;
ax2(m) = subplot(NRows,NCols,13,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(2,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1.5 1.5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'ytick',[-1.5 0 1.5],'yticklabel',{'-1.5','0','1.5'})
text(1.2, 0, '$x_2(n)$','interpreter','latex','fontsize',FS)

m=m+1;
ax2(m) = subplot(NRows,NCols,17,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(3,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1 1])
set(ax2(m),'fontsize',FS,'ytick',[-0.5 0 0.5],'yticklabel',{'-0.5','0','0.5'},'xtick',[t_start t_start+(t_end-t_start)/2 t_end])
text(1.2, 0, '$x_3(n)$','interpreter','latex','fontsize',FS)

% plot the second feature

t_start = 5.0;
t_end = t_start + WindowSize;
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'r','edgecolor','none','parent',ax2(1),'facealpha',0.4)

s_start = round(t_start*Fs);
s_end = round(t_end*Fs);

m=m+1;
ax2(m) = subplot(NRows,NCols,6,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),y(s_start:s_end),'k','parent',ax2(m))
hold(ax2(m),'on')
plot(t(s_start:s_end),sum(C(:,s_start:s_end),1),'r','parent',ax2(m))
hold(ax2(m),'off')
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'r','edgecolor','none','parent',ax2(m),'facealpha',0.4)

m=m+1;
ax2(m) = subplot(NRows,NCols,10,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(1,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,14,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(2,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1.5 1.5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,18,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(3,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1 1])
set(ax2(m),'fontsize',FS,'yticklabel',[],'xtick',[t_start t_start+(t_end-t_start)/2 t_end])

% plot the third feature

t_start = 8.5;
t_end = t_start + WindowSize;
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'b','edgecolor','none','parent',ax2(1),'facealpha',0.4)

s_start = round(t_start*Fs);
s_end = round(t_end*Fs);

m=m+1;
ax2(m) = subplot(NRows,NCols,7,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),y(s_start:s_end),'k','parent',ax2(m))
hold(ax2(m),'on')
plot(t(s_start:s_end),sum(C(:,s_start:s_end),1),'r','parent',ax2(m))
hold(ax2(m),'off')
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'b','edgecolor','none','parent',ax2(m),'facealpha',0.4)

m=m+1;
ax2(m) = subplot(NRows,NCols,11,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(1,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,15,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(2,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1.5 1.5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,19,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(3,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1 1])
set(ax2(m),'fontsize',FS,'yticklabel',[],'xtick',[t_start t_start+(t_end-t_start)/2 t_end])

% plot the forth feature
t_start = 12.1;
t_end = t_start + WindowSize;
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'y','edgecolor','none','parent',ax2(1),'facealpha',0.4)

s_start = round(t_start*Fs);
s_end = round(t_end*Fs);

m=m+1;
ax2(m) = subplot(NRows,NCols,8,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),y(s_start:s_end),'k','parent',ax2(m))
hold(ax2(m),'on')
plot(t(s_start:s_end),sum(C(:,s_start:s_end),1),'r','parent',ax2(m))
hold(ax2(m),'off')
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])
patch([t_start, t_end, t_end, t_start, t_start], [min(y), min(y), max(y), max(y), min(y)],'y','edgecolor','none','parent',ax2(m),'facealpha',0.4)

m=m+1;
ax2(m) = subplot(NRows,NCols,12,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(1,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-5 5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,16,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(2,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1.5 1.5])
set(ax2(m),'fontsize',FS,'xticklabel',[],'yticklabel',[])

m=m+1;
ax2(m) = subplot(NRows,NCols,20,'parent',TimeSeriesFig); 
plot(t(s_start:s_end),C(3,s_start:s_end),'k','parent',ax2(m))
xlim(ax2(m),[t_start t_end])
ylim(ax2(m),[-1 1])
set(ax2(m),'fontsize',FS,'yticklabel',[],'xtick',[t_start t_start+(t_end-t_start)/2 t_end])
