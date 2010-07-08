clc
clear
close all

Fs = 5e3;
Ts = 1/Fs;
T_max = 0.5;
t= 0:Ts:T_max-Ts;                   % 10 seconds @ 1kHz sample rate
f0 = 5;
f1 = 150;
a = linspace(2,0.5,length(t));
% y = a.*chirp(t,f0,T_max,f1) + 0.0*randn(1,length(t));
FS = 20;
LW = 1.5;
LW2 = 1.3;
ExtraWidth = 0.04;
Offset = 0.025;
ShiftUp = 0.05;

CuspFigure = figure('units','normalized','position',[0 0 1 0.4]);
a1 = [2 3 4];    
a2 = 1;
f1 = 20;
f2 = 60;   
g = abs((a1.^2*f1 + a2^2*f2)./(a1*a2*(f1+f2)));
N_imfs = 1;

for n=1:length(a1)
    
    y = a1(n)*cos(2*pi*20*t) + a2*cos(2*pi*60*t) + 0.0000*randn(1,length(t));
    Hy = a1(n)*sin(2*pi*20*t) + a2*sin(2*pi*60*t);
    ax_index = 1+3*(n-1);
    ax(ax_index) = subplot(length(a1),3,ax_index);
    plot(t,y,'k--','parent', ax(ax_index),'linewidth',LW)
    xlim(ax(ax_index),[0.1 0.4])
    ylim(ax(ax_index),[-6 6])
    hold(ax(ax_index),'on')
    
    InitialPoints = 4;
    PointsStep = 1;
    UpperLimit = 4;
    PlotMode = 0;

    % [IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0,
    % m_star, M] = cpt_rework(y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);
    [C r_approx IF_interp phi_interp phi_unwrapped m_star, M] = CPT_EMD_rework(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);
    
    plot(t,C,'r','parent',ax(ax_index),'linewidth',LW2)
%     plot(t,cos(phi_interp),'g','parent',ax(ax_index),'linewidth',LW2)
    
    if n<length(a1)
        set(ax(ax_index),'xticklabel',[],'fontsize',FS)
    else
        xlabel(ax(ax_index),'Time (s)','fontsize',FS)%,'interpreter','latex')
    end
    set(ax(ax_index),'fontsize',FS)
    ylabel(ax(ax_index),'$y(n)$','interpreter','latex','fontsize',FS)
    Pos = get(ax(ax_index),'position');
    if n==2
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2)+ShiftUp,Pos(3)+ExtraWidth,Pos(4)])
    elseif n>2
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2)+2*ShiftUp,Pos(3)+ExtraWidth,Pos(4)])
    else
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2),Pos(3)+ExtraWidth,Pos(4)])        
    end
    TextSpace = 2;
    text(-.03, TextSpace, ['$a_1 = ' num2str(a1(n)) '$'],'interpreter','latex','fontsize',FS)%, 'position',[0.5 0.5 0.1 0.1],'parent',CuspFigure)
    text(-.03, -TextSpace, '$a_2 = 1$','interpreter','latex','fontsize',FS)%, 'position',[0.5 0.5 0.1 0.1],'parent',CuspFigure)
    if n==1
        text(.222, 8, '$Q = 4$','interpreter','latex','fontsize',FS)%, 'position',[0.5 0.5 0.1 0.1],'parent',CuspFigure)
    end
    
    InitialPoints = 120;
    PointsStep = 1;
    UpperLimit = 120;
    PlotMode = 0;

    % [IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0, m_star, M] = cpt_rework(y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);

    [C r_approx IF_interp phi_interp phi_unwrapped m_star, M] = CPT_EMD_rework(N_imfs, y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode);
    
    ax_index = 2+3*(n-1);
    ax(ax_index) = subplot(length(a1),3,ax_index);
    plot(t,y,'k--','parent', ax(ax_index),'linewidth',LW)
    xlim(ax(ax_index),[[0.1 0.4]])
    ylim(ax(ax_index),[-6 6])
    hold(ax(ax_index),'on')
    
    plot(t,C,'r','parent',ax(ax_index),'linewidth',LW2 )
    set(ax(ax_index),'yticklabel',[],'fontsize',FS)
    if n<length(a1)
        set(ax(ax_index),'xticklabel',[],'fontsize',FS)
    else
        xlabel(ax(ax_index),'Time (s)','fontsize',FS)%,'interpreter','latex'
    end
    Pos = get(ax(ax_index),'position');
    if n==2
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2)+ShiftUp,Pos(3)+ExtraWidth,Pos(4)])
    elseif n>2
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2)+2*ShiftUp,Pos(3)+ExtraWidth,Pos(4)])
    else
        set(ax(ax_index),'position',[Pos(1)+Offset,Pos(2),Pos(3)+ExtraWidth,Pos(4)])
    end
    if n==1
        text(.211, 8, '$Q = 110$','interpreter','latex','fontsize',FS)%, 'position',[0.5 0.5 0.1 0.1],'parent',CuspFigure)
    end
    
    ax_index = n*3;
    ax(ax_index) = subplot(length(a1),3,ax_index);
    plot(y,Hy,'k','parent',ax(ax_index),'linewidth',LW2)
    xlim(ax(ax_index),[-6 6]), ylim(ax(ax_index),[-6 6])
    set(ax(ax_index),'fontsize',FS)
    axis(ax(ax_index),'square')
    set(ax(ax_index),'yticklabel',[],'fontsize',FS)
    if n<length(a1)
        set(ax(ax_index),'xticklabel',[],'fontsize',FS)
    else
        xlabel(ax(ax_index),'$y(n)$','interpreter','latex','fontsize',FS)
    end
    ylabel('$\mathcal{H} y(n)$','interpreter','latex','fontsize',FS)
    Pos = get(ax(ax_index),'position');
    if n==2
        set(ax(ax_index),'position',[Pos(1),Pos(2)+ShiftUp,Pos(3),Pos(4)])
    elseif n>2
        set(ax(ax_index),'position',[Pos(1),Pos(2)+2*ShiftUp,Pos(3),Pos(4)])
%     else
%         set(ax(ax_index),'position',[Pos(1),Pos(2)+ShiftUp,Pos(3),Pos(4)])
    end
end
% linkaxes(ax([1,4,7,10]),'x')
