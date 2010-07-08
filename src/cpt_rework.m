
function [IF_interp, phi_interp, phi_unwrapped, x, Hx, r, phi, x0, Hx0, m_star, M] = cpt_rework(y, InitialPoints, UpperLimit, PointsStep, Ts, PlotMode)

MS = 6;

Max_y = 2*max(y);
Min_y = 2*min(y);

XLimits = [Min_y Max_y];
YLimits = [Min_y Max_y];

Samples_For_Tangent = InitialPoints;       % must be an odd number

% we assume y is a NSamples x 1 vector
[NSamples NCols] = size(y);
if NCols > NSamples
    NSamples = NCols;
    y = y';
end

% AntiGibbsWindow = tukeywin(NSamples,0.05);          % use a tukey window to reduce edge effects
% y = AntiGibbsWindow.*y;                                            % make sure the signal is zero around the edges.

Hy = imag(hilbert(y));                          % compute the hilbert transform

x0 = nan*ones(1,NSamples);             % initialize translation parameters for speed
Hx0 = nan*ones(1,NSamples);            % there is a row for each direction
r = nan*ones(1,NSamples);              % there columns are the results
x = nan*ones(1,NSamples);
Hx = nan*ones(1,NSamples);
phi = nan*ones(1,NSamples);

e_star = zeros(1,NSamples);
d_star = zeros(1,NSamples);
m_star = zeros(1,NSamples);
flag = zeros(1,NSamples);
M = nan*ones(1,NSamples);    

StartSample = floor(InitialPoints/2) + 1;           % this is the first sample where we can get the number of points for the CPT
EndSample = NSamples-floor(InitialPoints/2);    % this is the last sample the we can get the CPT given the smallest arc length

flag(1:StartSample) = 1;                                    % to indicate this are not good mappings

for n=StartSample:EndSample
    
    if n < NSamples/2               % This is not what is used, but what is available. Work out the largest number of points available to create mapping
        MaxPoints = n*2-1;          % at the beginning and at the end of the time series there is fewer points available
    else
        MaxPoints = 2*(NSamples-n)+1;
    end

    if MaxPoints > UpperLimit
        MaxPoints = UpperLimit;     % we are far enough away from the ends of the time series to use upper limit on number of points to create map
    end    
    
    PointsFromCenter = floor(MaxPoints/2);
    
    y_s_all = y(n-PointsFromCenter:n+PointsFromCenter);                 % take our segment of data
    Hy_s_all = Hy(n-PointsFromCenter:n+PointsFromCenter);               % the seg size is the maximum for as specified for this data point
    
    [z a] = Get_Tangent_and_Normal_Vectors(Samples_For_Tangent,y_s_all, Hy_s_all);   
    
    if PlotMode == 1
        plot(y_s_all,Hy_s_all,'k*')
        hold on
        plot(y(n),Hy(n),'r*')
%         plot([0 z(1)],[0 z(2)],'r')
%          plot([0 a(1)],[0 a(2)],'k')
        xlim(XLimits)
        ylim(YLimits)
        axis square
    end
    
    x_previous = x(n-PointsFromCenter:n-1);                                         % note: current sample not include! take previous estimates corresponding to segment size
    Hx_previous = Hx(n-PointsFromCenter:n-1);                                     % this is used for the relative distance mismatch
      
    [x0_temp Hx0_temp r_temp phi_temp x_temp Hx_temp M_temp e d Norm] =  Find_All_Circle_Parameters_And_Tangents(MaxPoints, InitialPoints, PointsStep, y_s_all, Hy_s_all, x_previous, Hx_previous);
    
    if PlotMode == 1
        hold on
        plot(x0_temp,Hx0_temp,'+')
        hold off
        xlim(XLimits)
        ylim(YLimits) 
        axis square
    end
    
    % if lambda < 0 then circle is on the correct side
    lambda = Check_Circle_Center(y(n), Hy(n), x0_temp, Hx0_temp, x_temp, Hx_temp, a);
%     lambda = Check_Circle_Center(y(n), Hy(n), x0_temp, Hx0_temp, x_temp, Hx_temp, Norm);
    
    if sum(lambda<0) == 0
        flag(n) = 1;                    % means no circles were on the correct side
        if PlotMode == 2 || PlotMode == 1
            disp('All circle fits were on the incorrect side of the estimated tangent')
        end
    else
        
        % restrict mapping that have have the centers on the correct side of
        % the estimated tangent    
        x0_temp = x0_temp(lambda<0);
        Hx0_temp = Hx0_temp(lambda<0); 
        r_temp  = r_temp(lambda<0);
        phi_temp = phi_temp(lambda<0);
        x_temp = x_temp(lambda<0); 
        Hx_temp = Hx_temp(lambda<0); 
        M_temp = M_temp(lambda<0); 
        e = e(lambda<0); 
        d = d(lambda<0);

        last_good_mapping_index = find(flag(1:n-1)==0,1,'last');        % we will use this index with our priors
        
        if isempty(last_good_mapping_index)
            
            [e_star(n) m_star(n)] = min(abs(e));
            x0(n) = x0_temp(m_star(n));      % these are the best guess of the circle parameters
            Hx0(n) = Hx0_temp(m_star(n));
            r(n) = r_temp(m_star(n));        
            phi(n) = phi_temp(m_star(n));
            x(n) = x_temp(m_star(n));
            Hx(n) = Hx_temp(m_star(n));
            M(n) = M_temp(m_star(n));
            
        else
        
            [p1 p2] = Calculate_Priors(r_temp, phi_temp, r(last_good_mapping_index), phi(last_good_mapping_index));

            PhaseForwardThresh = 2*pi;%.225*pi;        %0.17
            PhaseJumpThresh = -0.75*pi;        % this is used with the condtion that the phase difference phi(n)-phi(n-1) is positive. It can also be less than this consant to account for 2pi phase jumps.

%             Indexes_From_p1 = ((p1<PhaseForwardThresh) & (p1 > 0)) | (p1<PhaseJumpThresh);
            Indexes_From_p1 = ((p1<PhaseForwardThresh) & (p1 >= 0)) | (p1<PhaseJumpThresh);
            Indexes_From_p1_and_p2 = (((p1<PhaseForwardThresh) & (p1 >= 0)) | (p1<PhaseJumpThresh)) & (p2>0);

            if sum(Indexes_From_p1) > 0         % if this is not satisfied than we have to assign NaN

                if sum(Indexes_From_p1_and_p2) > 0     % if this is satisfied than we will use it, otherwise we will not

                    if PlotMode == 2 
                        disp('p1 satisfied, p2 satisfied')
                    end         

                    x0_temp = x0_temp(Indexes_From_p1_and_p2 );
                    Hx0_temp = Hx0_temp( Indexes_From_p1_and_p2); 
                    r_temp  = r_temp( Indexes_From_p1_and_p2);
                    phi_temp = phi_temp(Indexes_From_p1_and_p2 );
                    x_temp = x_temp(Indexes_From_p1_and_p2 ); 
                    Hx_temp = Hx_temp( Indexes_From_p1_and_p2); 
                    M_temp = M_temp( Indexes_From_p1_and_p2 ); 
                    e = e( Indexes_From_p1_and_p2 ); 
                    d = d( Indexes_From_p1_and_p2 );

                else

                    if PlotMode == 2 
                        disp('p1 satisfied, p2 not satisfied')
                    end

                    x0_temp = x0_temp( Indexes_From_p1);
                    Hx0_temp = Hx0_temp(Indexes_From_p1 ); 
                    r_temp  = r_temp( Indexes_From_p1);
                    phi_temp = phi_temp(Indexes_From_p1 );
                    x_temp = x_temp(Indexes_From_p1 ); 
                    Hx_temp = Hx_temp( Indexes_From_p1); 
                    M_temp = M_temp( Indexes_From_p1 ); 
                    e = e(Indexes_From_p1  ); 
                    d = d( Indexes_From_p1 );

                end

                if sum(~isnan(d))==0            % then we cant use d and we have to use e
                    if PlotMode == 2
                        disp('d is empty, using e')
                    end
                    [e_star(n) m_star(n)] = min(abs(e));
                else
                    [d_star(n) m_star(n)] = min(abs(d));
                end  

                x0(n) = x0_temp(m_star(n));      % these are the best guess of the circle parameters
                Hx0(n) = Hx0_temp(m_star(n));
                r(n) = r_temp(m_star(n));        
                phi(n) = phi_temp(m_star(n));
                x(n) = x_temp(m_star(n));
                Hx(n) = Hx_temp(m_star(n));
                M(n) = M_temp(m_star(n));

            else

                if PlotMode == 2 
                    disp('p1 not satisfied')
                end
                flag(n) = 1;

            end
        end
    end
    
    if PlotMode == 1
        if ~isnan(M(n))
            hold on
            plot(x(n-floor(M(n)/2):n),Hx(n-floor(M(n)/2):n),'go','markersize',MS)
            hold off
            xlim(XLimits)
            ylim(YLimits)
            drawnow
            axis square
%         pause(0.2)
        else
%             disp('M(n) is nan')
        end
    end
    
end

% unwrap phase
phi_unwrapped = unwrap_phi(phi);

% interpolate missing values
[IF_interp phi_interp] = interpolate_IP_and_IF(phi_unwrapped,Ts);
