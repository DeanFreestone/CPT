
% emd cpt function 
% Dean Freestone 13/01/2010
% UpperLimit vector same length as y all the same eg.
% UpperLimit=10*ones(1,length(y));

% C = IMFs
% r_approx = approximate IA (emperical)
% phi_interp = IP with NaNs interpolated
% phi_unwrapped = with NaNs

function [C r_approx IF_interp phi_interp phi_unwrapped m_star, M, res,x0,Hx0] = ...
    CPT_EMD(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode)

[NCols NSamples] = size(y);
if NSamples < NCols
    y = y';
end

NSamples = length(y);
phi_interp = zeros(N_imfs,NSamples);
phi_unwrapped = zeros(N_imfs,NSamples);
r_approx = zeros(N_imfs,NSamples);
C = zeros(N_imfs,NSamples);
res = zeros(N_imfs,NSamples);
M = zeros(N_imfs,NSamples);
s1 = zeros(N_imfs,NSamples);
s2 = zeros(N_imfs,NSamples);

yy = y;

for n=1:N_imfs
    
    disp(['Processing ' num2str(n) ' of ' num2str(N_imfs) ' IMFs'])
    
    if n>1
        InitialPoints = InitialPoints;
        UpperLimit = UpperLimit;
    end
        
    [IF_interp(n,:), phi_interp(n,:), phi_unwrapped(n,:), x, Hx, r, phi, x0, Hx0, m_star(n,:), M(n,:)] = cpt(yy,InitialPoints,UpperLimit, PointsStep, Ts, PlotMode,f_max);
    
    % find the detrending line
    cos_phi = cos(phi_interp(n,:));
    p1 = findpeaks(cos_phi);
    p2 = findpeaks(-cos_phi);
    p3 = sort([p1 p2]);
    
    % still not sure if pchip is better!
    s1(n,:) = spline([0 p1 NSamples+1],[0 yy(p1) 0],1:NSamples);                  % need to check if pchip works better
    s2(n,:) = -spline([0 p2 NSamples+1],[0 -yy(p2) 0],1:NSamples);                % also need to check if 1 and NSamples works better for indexes

    res(n,:) = (s1(n,:)+s2(n,:))/2;        % this is the residue or the detrending line
    
%     sin_phi = sin(phi_interp(n,:));
%     p11 = findpeaks(sin_phi);
%     p22 = findpeaks(-sin_phi);
%     
%     % still not sure if pchip is better!
%     s11(n,:) = spline([0 p11 NSamples+1],[0 yy(p11) 0],1:NSamples);                  % need to check if pchip works better
%     s22(n,:) = -spline([0 p22 NSamples+1],[0 -yy(p22) 0],1:NSamples);                % also need to check if 1 and NSamples works better for indexes
% 
%     res2(n,:) = (s11(n,:)+s22(n,:))/2;        % this is the residue or the detrending line
%     
    
    % detrend and find envelope
    detrend = yy-res(n,:);
    abs_detrend = abs(detrend);
    r_approx(n,:) = spline([0 p3 NSamples+1],[0 abs_detrend(p3) 0],1:NSamples);   % this is an approximate envelope

    % assign IMF
    C(n,:) = r_approx(n,:).*cos(phi_interp(n,:));

    plotmode =0;
    t = linspace(0,3,length(yy));

    if (n == 1) && (plotmode == 1)
       
        CPTEMDExample = figure('name','decomposition example','units','normalized','position',[0 0 0.45 0.6],...
            'PaperOrientation', 'landscape','papertype','a4','paperunits','normalized');%,'PaperPositionMode','auto')
        
        FS = 16;
        FSleg = 15;
        MS = 16;
        LW1 = 3;
        ylims = [-.3 .3];
        tlims = [1 2];
        
        axe(1) = subplot(311);
        plot(t,yy, 'linewidth', LW1)
        hold
        
        plot(t,s1(n,:),'-.g', 'linewidth', LW1)
        plot(t,s2(n,:),'-.c', 'linewidth', LW1)
        plot(t,res(n,:),'-.k', 'linewidth', LW1)
        plot(t(p3),yy(p3),'.r','markersize',MS)
        hold off
        set(gca,'xticklabel',[],'fontsize',FS)
        ylabel('Amplitude','fontsize',FS)%,'interpreter','latex')
        title('A.','fontsize',FS,'units','normalized','position',[-0.1 1.1 0])% , 'interpreter','latex'       
        leg1 = legend('$y(n)$', '$s_{max,k }(n)$', '$s_{min,k }(n)$', '$x_{0,k}(n)$', 'extrema');
        set(leg1, 'interpreter', 'latex', 'fontsize',FSleg, 'orientation', 'horizontal', 'edgecolor', [1 1 1],...
            'position', [0.3 0.94 0.35 0.030])   
        ylim(ylims)
        xlim(tlims)
        
        axe(2) = subplot(312);
        plot(t,abs_detrend, 'linewidth', LW1)
        hold
        plot(t,r_approx(n,:),'-.g', 'linewidth', LW1)
        plot(t(p3),abs_detrend(p3),'.r','markersize',MS)
        hold off
        set(gca,'xticklabel',[],'fontsize',FS)
        ylabel('Amplitude','fontsize',FS)%,'interpreter','latex')
        leg2 = legend('$h_k(n)$', '$r_k(n)$', 'maxima');
        set(leg2, 'interpreter', 'latex', 'fontsize',FSleg, 'orientation', 'horizontal', 'edgecolor', [1 1 1],...
            'position', [0.27 0.645 0.5 0.04])
        title('B.','fontsize',FS, 'units','normalized','position',[-0.1 1.1 0])%,'interpreter','latex')
        ylim([0 ylims(2)])
        xlim(tlims)
        
        axe(3) = subplot(313);
        plot(t,r_approx(n,:).*cos_phi, 'linewidth', LW1)
        xlabel('Time (s)','fontsize',FS)%,'interpreter','latex')
        ylabel('Amplitude','fontsize',FS)%,'interpreter','latex')
        title('C.','fontsize',FS,'units','normalized','position',[-0.1 1.1 0])%, 'interpreter','latex'
        set(gca,'fontsize', FS)
        ylim(ylims)
        leg3 = legend('$x_k(n)$');
        set(leg3, 'interpreter', 'latex', 'fontsize',FSleg, 'orientation', 'horizontal', 'edgecolor', [1 1 1],...
            'position', [0.44 0.34 0.155 0.04])
        xlim(tlims)
        
        linkaxes(axe,'x')
        
        % Set the page size and position to match the figure's dimensions
        set(CPTEMDExample,'PaperUnits','inches');
        set(CPTEMDExample,'Units','inches');
        paperPosition = get(CPTEMDExample,'PaperPosition');
        position = get(CPTEMDExample,'Position');
        set(CPTEMDExample,'PaperPosition',[0,0,position(3:4)]);
        set(CPTEMDExample,'PaperSize',position(3:4));
        
%         print(CPTEMDExample,'-dpdf','CPTEMDExample.pdf',sprintf('-r%d',600))
    elseif plotmode == 2
        
        PlotFunction151009
        
    end
    
    % for next iteration
    yy=res(n,:);
end