clear
clc
close all

% FS = 18;                % fontsize for figures
FS2 = 10;               % for the legend
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

MinPoints = round(0.7*Fs/f1);
MaxPoints = round(PointMultiplier*MinPoints);
PointsStep = 2;
PlotMode = 0;

PhaseFig  = figure('units','normalized','position',[0 0 1 1]);
FS = 16;
LW = 1.6;
MS = 5;

EdgeEffectTime = 0.02;
EdgeEffectSample = round(EdgeEffectTime*Fs);
t_start_sample = 1;
NIterations = 100;

RMSE_start = zeros(NIterations,length(sigma));
RMSE_mid = zeros(NIterations,length(sigma));
RMSE_end = zeros(NIterations,length(sigma));

tic
for n=1:length(sigma)
    disp(['iteration ' num2str(n) ' of ' num2str(length(sigma)) ' noise levels']) 
    
    for nn=1:NIterations
        
        noise = sigma(n)*randn(1,length(t));
        y_n = y + noise;

        s_1 = round(Fs*T_max/3);
        s_2 = round(2*Fs*T_max/3);

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

        phi_est = mod(phi_interp,2*pi);

        phi_diff = min( (repmat(true_IP,3,1) - [phi_est-2*pi ; phi_est ; phi_est + 2*pi]).^2 ,[],1); 
        phi_diff_start = phi_diff(EdgeEffectSample:s_1);
        phi_diff_mid = phi_diff(s_1:s_2);
        phi_diff_end = phi_diff(s_2:end-EdgeEffectSample);

        RMSE_start(nn,n) = sqrt(mean(phi_diff_start));
        RMSE_mid(nn,n) = sqrt(mean(phi_diff_mid));
        RMSE_end(nn,n) = sqrt(mean(phi_diff_end));


    %     phi_diff_mid
    %     phi_diff_end
% 
%         ax(n+floor((n-1)/3)*3) = subplot(6,3,n+floor((n-1)/3)*3,'parent',PhaseFig);
%         plot(t,y_n,'k','parent',ax(n+floor((n-1)/3)*3))
%         grid(ax(n+floor((n-1)/3)*3),'on')
%         Pos = get(ax(n+floor((n-1)/3)*3),'position');
%         set(ax(n+floor((n-1)/3)*3),'position',[Pos(1) Pos(2) Pos(3) Pos(4)-0.01]);
%         title(ax(n+floor((n-1)/3)*3),['SNR = ' num2str(snr_db_start) ', ' num2str(snr_db_mid) ', ' num2str(snr_db_end) ' dB'],'fontsize',FS,'fontname','Times New Roman')%'interpreter','latex')
%         if mod(n,3) ~= 1
%             set(ax(n+floor((n-1)/3)*3),'yticklabel',[])
%         else
%             ylabel(ax(n+floor((n-1)/3)*3),'Amplitude','fontsize',FS,'fontname','Times New Roman','position',[-0.031 -0.5 1])%,'interpreter','latex')
%         end
%         set(ax(n+floor((n-1)/3)*3),'xticklabel',[],'fontsize',FS,'fontname','Times New Roman','ytick',[-2 0 2 4])
%         xlim(ax(n+floor((n-1)/3)*3),[t(t_start_sample) t(end)])
%         ylim(ax(n+floor((n-1)/3)*3),[-A1 A1])
% 
%         ax(6,3,n+(floor((n-1)/3)+1)) = subplot(6,3,n+(floor((n-1)/3)+1)*3,'parent',PhaseFig);
%         plot(t,true_IP,'k','linewidth',LW,'parent',ax(6,3,n+(floor((n-1)/3)+1)))
%         hold(ax(6,3,n+(floor((n-1)/3)+1)),'on')
%         plot(t,phi_est,'.r','markersize',MS,'parent',ax(6,3,n+(floor((n-1)/3)+1)))
%         hold(ax(6,3,n+(floor((n-1)/3)+1)),'off')
%         Pos = get(ax(6,3,n+(floor((n-1)/3)+1)),'position');
%         set(ax(6,3,n+(floor((n-1)/3)+1)),'position',[Pos(1) Pos(2)+0.01 Pos(3) Pos(4)+0.02]);
%         grid(ax(6,3,n+(floor((n-1)/3)+1)),'on')
%         if n<7
%             set(ax(6,3,n+(floor((n-1)/3)+1)),'xticklabel',[])
%         else
%             xlabel(ax(6,3,n+(floor((n-1)/3)+1)),'Times (s)','fontsize',FS,'fontname','Times New Roman')%,'interpreter','latex')
%         end
%         set(ax(6,3,n+(floor((n-1)/3)+1)),'fontsize',FS,'fontname','Times New Roman','ytick',[0 pi 2*pi],'yticklabel',{'0', '3.14','6.28'})
%         if mod(n,3) ~= 1
%             set(ax(6,3,n+(floor((n-1)/3)+1)),'yticklabel',[])
%         else
%             ylabel(ax(6,3,n+(floor((n-1)/3)+1)),'Radians','fontsize',FS,'fontname','Times New Roman','position',[-0.031 3 1])%,'interpreter','latex')
%         end
%         xlim(ax(6,3,n+(floor((n-1)/3)+1)),[t(t_start_sample) t(end)])
%         ylim(ax(6,3,n+(floor((n-1)/3)+1)),[0,2*pi+0.4])
%         drawnow

    end
end
toc

for n=1:length(sigma)
    disp(['iteration ' num2str(n) ' of ' num2str(length(sigma)) ' noise levels']) 
    
    for nn=1:NIterations
        
        noise = sigma(n)*randn(1,length(t));
        y_n = y + noise;

        s_1 = round(Fs*T_max/3);
        s_2 = round(2*Fs*T_max/3);

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
    end
end

y_max= 1.5;
FS = 16;
roundlabelstart = round(10*mean(snr_db_start,1))/10;
roundlabelmid = round(10*mean(snr_db_mid,1))/10;
roundlabelend = round(10*mean(snr_db_end,1))/10;
for n=1:length(sigma)
    xticklabelstart{n} = num2str(roundlabelstart(n));
    xticklabelmid{n} = num2str(roundlabelmid(n));
    xticklabelend{n} = num2str(roundlabelend(n));
end
boxplotfig = figure('units','normalized','position',[0 0 1 1]);
ax(1) = subplot(3,3,1,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_start,roundlabelstart,'parent',ax(1))
ylim(ax(1),[0 y_max])
pos = get(ax(1),'position');
set(ax(1),'xtick',1:9,'xticklabel',xticklabelstart,'fontsize',FS,'position',[pos(1) pos(2) pos(3)+0.01 pos(4)])

ax(2) = subplot(3,3,2,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_mid,roundlabelmid,'parent',ax(2))
ylim(ax(2),[0 y_max])
pos = get(ax(2),'position');
set(ax(2),'xtick',1:9,'xticklabel',xticklabelmid,'fontsize',FS,'position',[pos(1) pos(2) pos(3)+0.01 pos(4)])

ax(3) = subplot(3,3,3,'parent',boxplotfig,'units','normalized');
boxplot(RMSE_end,roundlabelend,'parent',ax(3))
ylim(ax(3),[0 y_max])
pos = get(ax(3),'position');
set(ax(3),'xtick',1:9,'xticklabel',xticklabelend,'fontsize',FS,'position',[pos(1) pos(2) pos(3)+0.01 pos(4)])

% clims_db = [-3 27];
% clims_rmse = [0 0.9];
% 
% figure('units','normalized','position',[0 0 1 1])
% subplot(231)
% imagesc(PointMultiplier,sigma,mean(snr_db,3),clims_db)
% axis square
% colorbar
% title('SNR db Original','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% subplot(234)
% imagesc(PointMultiplier,sigma,mean(rmse,3),clims_rmse)
% axis square
% title('RMSE Original','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% colorbar
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% subplot(232)
% imagesc(PointMultiplier,sigma,mean(snr_db_cpt,3),clims_db)
% axis square
% colorbar
% title('SNR db CPT','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% subplot(235)
% imagesc(PointMultiplier,sigma,mean(rmse_cpt,3),clims_rmse)
% axis square
% title('RMSE CPT','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% colorbar
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% subplot(233)
% imagesc(PointMultiplier,sigma,mean(rmse_cpt_IP,3))
% axis square
% title('RMSE IP CPT','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% colorbar
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% subplot(236)
% imagesc(PointMultiplier,sigma,mean(rmse_cpt_IF,3))
% axis square
% title('RMSE IF CPT','fontsize',FS)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('$\sigma$','fontsize',FS3,'interpreter','latex')
% colorbar
% set(gca,'fontsize',FS,'xtick',PointMultiplier(1:2:end),'ytick',sigma(1:2:end))
% 
% diff_snr = mean(snr_db,3)-mean(snr_db_cpt,3);
% figure
% imagesc(PointMultiplier,sigma,diff_snr)

% % ~~~~~~~~~~~
% figure('units','normalized','position',[0 0 1 1])
% ax3(1) = subplot(211);
% plot(PointMultiplier,squeeze(mean(snr_db_cpt,3)),'-*')
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% ylabel('SNR dB','fontsize',FS3)
% for n=1:length(sigma)
%     leg_entry{n} = ['$\sigma = ' num2str(sigma(n)) '$'];
% end
% leg1 = legend(leg_entry);
% set(leg1,'orientation' ,'horizontal','fontsize',FS2,'Position',[0.0407986111111113 0.951466916354557 0.931944444444444 0.0282459425717853],'interpreter','latex');
% xlim([PointMultiplier(1) PointMultiplier(end)])
% set(gca,'fontsize',FS,'xtick',PointMultiplier)
% 
% ax3(2) = subplot(212);
% plot(PointMultiplier,squeeze(mean(rmse_cpt,3)),'-*')
% ylabel('RMSE','fontsize',FS3)
% xlabel('$\alpha$','fontsize',FS3,'interpreter','latex')
% set(gca,'fontsize',FS,'xtick',PointMultiplier)
% xlim([PointMultiplier(1) PointMultiplier(end)])