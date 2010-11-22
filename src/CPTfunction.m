function [foundarc phi phi_unwrapped r t firstindex lastindex ArcPoints TangentPoints] = CPTfunction(y, t, Ts, psi, f_max,init_b_f,zeta)

% ploton = true;
ploton = false;

CounterLimit = 5; 

% initialize things for speed
foundarc = false;                       % this is a flag that is used in the cpt emd. it is set if an arc is found. if not the emd should end.
NSamples = length(y);
TangentPoints = zeros(2,NSamples);
ArcPoints = zeros(2,NSamples);
Tangent = zeros(2,NSamples);
a = zeros(2,NSamples);
x0 = zeros(1,NSamples);                                     % initialize estimate for x0 translation
Hx0 = zeros(1,NSamples);                                    % initialize estimate for Hx0 translation
r_temp = zeros(1,NSamples);                                % initialize estimate for radius (envelope)
r = zeros(1,NSamples);
phi = zeros(1,NSamples);

Hy = imag(hilbert(y));                                  % Hilbert transform of y
P_y = [y ; Hy];
                                          % this counter is for detecting cusps and the desired arc length
SampleThreshold = 2/Ts;

SampleForwardBackFlag = false;              % if the number of sample to fit arc gets too big, change flag

init_b = init_b_f;     
init_f = init_b_f;      

% cycle through all samples to find circle fits
for n=init_b(1)+1:NSamples - init_f(1)-zeta
    
    if ~SampleForwardBackFlag               % this means we have had to look too far to find an arc, so we stop. If we dont reset the flag we dont go into the while loop.
        flag = false;                                                                   % initialize / reset arc length indicator flag, this flag tells us when to get out of the while loop below
    end
    
    P_yn = P_y(:,n);                                                            % point of interest where we want to fit the circle
    
    % now there is no guarantee for the emd that the n-b or n+f will be in
    % the data range because the signal gets shorter with each step of the
    % decomposition. Therefore we need to make some checks
    
    % first check if we are in range
    if (init_b(n) < n) && (n+init_f(n)  < NSamples)
        b_t_final = init_b(n);                         % number of samples back in time from point of interest to build arc for circle fit
        f_t_final = init_f(n);                          % number of samples forward in time
        
     % if we are not within range check if we are over the edge and make
     % the size of the arc as large as possible
    elseif (n+init_f(n)  >= NSamples)
        b_t_final = floor((NSamples - n)/2)-1;
        f_t_final = b_t_final;
        flag = true;                    % we don't need to go through the loop because we at an edge
    else
        
    % we are over the start, so we make the arc as large as possible.
        b_t_final = n-1;                         % number of samples back in time from point of interest to build arc for circle fit
        f_t_final = n-1;                            % number of samples forward in time
        flag = true;                    % we don't need to go through the loop because we at an edge
    end
    
    b = b_t_final;
    f = f_t_final;
    b_final = b;
    f_final = f;

%     theta = pi;                                     % initialize angle describing arc length, which should be decreasing as we iterate through the while loop
    ArcCounter = 0;                             % this counter helps us deal with noise by making sure we have a suficent number of arc below the desired arc length to be sure noise hasn't caused a false detection
    TangentCounter = 0;
    CuspCounter = 0;
    
    % we go through this only on the first time through
    p_b = P_y(:,n-b) - P_yn;
    p_f = P_y(:,n+f) - P_yn;
    
    theta = acos(dot(p_f,p_b) / (norm(p_b)*norm(p_f)));
    psi_hat = 2*(pi - theta);                               % estimated arc length
            
    % cycle through various arc lengths until the arc angle is correct
    while flag ~= true
               
        if (b > n) || (n+f > NSamples)
            disp('error')                           % for checking
        end
        
        % these distances will decide if we increment b or f.
        DistBack = norm(P_yn - P_y(:,n-b));            % Euclidean distance between center and back edge of arc
        DistForward = norm(P_yn - P_y(:,n+f));      % Euclidean distance between center and forward edge of arc
        
        if DistBack <= DistForward

            if n-round(b+zeta) > 0                         % make sure we don't go into negative indexes

                b = round(b + zeta); 
                p_b = P_y(:,n-b) - P_yn;

            else                                         % we are at the edge of the time series so we need to stop

                flag = true;                                        % we can't index any further back

                b_t_final = b;
                f_t_final = f;

                b_final = b;
                f_final = f;

            end

        elseif DistBack > DistForward

            if n+round(f+zeta) < NSamples           % make sure we don't try to index something longer than the time series

                f = round(f + zeta);
                p_f = P_y(:,n+f) - P_yn;

            else

                flag = true;                                        % we can't index any further forward

                b_t_final = b;
                f_t_final = f;

                b_final = b;
                f_final = f;

            end
        end
            
        % find angle between points at the edge of the arc
        theta_previous = theta;
        theta = acos(dot(p_f,p_b) / (norm(p_b)*norm(p_f)));
        psi_hat = 2*(pi - theta);                               % estimated arc length
                
        % check for cusps
        % ~~~~~~~~~
        if theta > theta_previous
            % means angle is opening up due to a cusp or noise
            if CuspCounter == 0
                b_final = b;
                f_final = f;
                b_t_final = b;
                f_t_final = f;
            end
            
            CuspCounter = CuspCounter + 1;
            
            if CuspCounter >= CounterLimit
                flag = true;                            % means we are sure we have found a cusp and we need to stop
            end
            
        else
            CuspCounter = 0;                    % reset cusp counter
        end
        
        % check for tangent
        % ~~~~~~~~~~~
        if psi_hat >= pi
            
            % then it might be the right points for the tangent
            
            if TangentCounter == 0
                b_t_final = b;
                f_t_final = f;
            end
            
            TangentCounter = TangentCounter + 1;
            
            if (TangentCounter >= CounterLimit)
                
                % then we are confident that it is the correct tangent
                % points, and not detected by noise
                
                if pi >= psi
                    
                    flag = true;                % this means we can stop because the desired arc length is less than the tangent
                    foundarc = true;        % means we dont need to stop and EMD

                end                
            end
        else
            
            % than we start looking again
            TangentCounter = 0;
            
        end

        % here we check if the arc length
        % ~~~~~~~~~~~~~~~~~~~
        if (psi_hat > psi) % || (theta_diff > theta_diff_thresh)
            
            if ArcCounter == 0
                % save then number of samples forward and back and what was
                % the increment...
                b_final = b;
                f_final = f;

            end
            
            ArcCounter = ArcCounter + 1;
            
            if ArcCounter >= CounterLimit
                
                % than this is our arc
                if psi >= pi
                    % than we don't need to find tangent
                    flag = true;        % we have the appropriate angle consistantly over 5 samples, therefore it is unlikely to be caused by noise
                end
                
                foundarc = true;            % means we can keep going with the EMD

            end
            
        else
            
            ArcCounter = 0;
            
        end
                
        % if b or f have gotten to big we want to stop.
        if (b > SampleThreshold) || (f > SampleThreshold)           
            SampleForwardBackFlag = true;
            flag = true;
        end
        
    end
    
    TangentPoints(:,n) = [b_t_final ; f_t_final];
    ArcPoints(:,n) = [b_final ; f_final];

    y_s = y( n-ArcPoints(1,n) : n+ArcPoints(2,n) );
    Hy_s = Hy( n-ArcPoints(1,n) : n+ArcPoints(2,n) );

    if n+TangentPoints(2,n) > NSamples
        disp('error')
    end
    Tangent(:,n) = [y(n+TangentPoints(2,n)) - y(n-TangentPoints(1,n)) 
        Hy(n+TangentPoints(2,n)) - Hy(n-TangentPoints(1,n))];

    a(:,n) = [Tangent(2,n) ; -Tangent(1,n)];      % normal to the tangent

    XY = [y_s ; Hy_s]';                                                % segment in matrix so we can use the circle fit

    Par = CircleFitByTaubin(XY);                    % Nonlinear. from Mathworks online. Uses LMS fit
%     Par = CircleFitByPratt(XY);                    % Nonlinear. from Mathworks online. Uses LMS fit

    x0(n) = Par(1);                                % estimate for y0 translation
    Hx0(n) = Par(2);                               % estimate for Hy0 translation
    r_temp(n) = Par(3);                                 % estimate for radius (envelope)

end

if ~SampleForwardBackFlag
    phi_temp = atan2(Hy-Hx0, y-x0);
    x_temp = r_temp.*cos(phi_temp);
    Hx_temp = r_temp.*sin(phi_temp);

    % need to check that the parameters are ok.
    % ~~~~~~~~~~~~~~~~~~~~~~~~

    % first check if the circle fit is on the correct side of the arc.
    lambda = zeros(1,NSamples);
    P_x = [x_temp; Hx_temp];
    for n=1:NSamples
        lambda(n) = a(:,n)'*P_x(:,n);
    end

    y_max = max(abs(y)-mean(y));

    % need to reject the bad estimates here before doing the phase bit
    % because I don' t want to include bad estimates
    CircleCenterIndexes = lambda<=0;
    YmaxIndexes = y_max < abs(x_temp);
    
    IndexesForNan = CircleCenterIndexes | YmaxIndexes;
    
    phi = phi_temp+pi;  % add pi to phase so it is between [0 2pi] so our checks are easier 
    
    r = r_temp;
    x = x_temp;
    Hx = Hx_temp;
    x(IndexesForNan) = nan;
    Hx(IndexesForNan) = nan;
    phi(IndexesForNan) = nan;
    r(IndexesForNan) = nan;
    
%     ~~~~~~~~~~~~~~~~~
    rho = 1/(Ts*2*f_max);           % this is the over sampling parameters
    p = 1;    
    beta_p = p*pi / rho;
    alpha_p = 2*beta_p - 2*pi;
        
    delta_phi = [0 phi(2:end) - phi(1:end-1)];
    
    % detect the invalid phase changes, find indexes for incorrect phase
    % transition
    PhaseIncIndexes = ~(( (delta_phi > 0) & (delta_phi < 2*beta_p) ) | (delta_phi < alpha_p));

    % we don't know if it was the later or earlier phase estimate that was
    % know so we set them both to nan.
%     temp = [PhaseIncIndexes(2:end) false];
%     PhaseIncIndexes = PhaseIncIndexes | temp;
    
    phi(PhaseIncIndexes) = nan; 
    x(PhaseIncIndexes) = nan;
    Hx(PhaseIncIndexes) = nan;
    r(PhaseIncIndexes) = nan;
%     ~~~~~~~


    % find the first IP estimate that is not nan
    first_good_index = find(~isnan(phi),1,'first');
    phi_previous = phi(first_good_index);
    rho = 1/(Ts*2*f_max);           % this is the over sampling parameters
    p = 1;
    for n=first_good_index+1:length(phi)
        phi_current = phi(n);
        if isnan(phi_current)
            p=p+1;
        else
            beta_p = p*pi / rho;
            alpha_p = 2*beta_p - 2*pi;
            delta_phi = phi_current - phi_previous;
            if ~(( (delta_phi > 0) && (delta_phi < 2*beta_p) ) || (delta_phi < alpha_p))
                phi(n) = nan;
                r(n) = nan;
                x(n) = nan;
                Hx(n) = nan;
                p = p + 1;
            else
                phi_previous = phi_current;
                p = 1;
            end
        end
    end

    if ploton
        figure
        PlotLim = 1.5*max(y);
        MS = 20;
        for n=b_init(1)+1:NSamples-f_init(1)-1

            y_s = y(n-ArcPoints(1,n):n+ArcPoints(2,n));
            Hy_s = Hy(n-ArcPoints(1,n):n+ArcPoints(2,n));
            plot(y_s,Hy_s,'ko-')
            hold on
            x_s_temp = x_temp(n-ArcPoints(1,n):n+ArcPoints(2,n));
            Hx_s_temp = Hx_temp(n-ArcPoints(1,n):n+ArcPoints(2,n));
            plot(x_s_temp,Hx_s_temp,'o-')

            x_s = x(n-ArcPoints(1,n):n+ArcPoints(2,n));
            Hx_s = Hx(n-ArcPoints(1,n):n+ArcPoints(2,n));
            plot(x_s,Hx_s,'g.-')

            plot(x0(n),Hx0(n),'*')

            plot(y(n),Hy(n),'xr','markersize',MS)

            plot([0 Tangent(1,n)],[0 Tangent(2,n)])
            plot([0 a(1,n)],[0 a(2,n)],'r')

            hold off
            xlim([-PlotLim PlotLim])
            ylim([-PlotLim PlotLim])
%             pause(0.05)
            drawnow

        end
    end
    % unwrap phase
    phi_unwrapped = unwrap_phi(phi) - pi;

    % interpolate missing values
    [IF_interp phi r t firstindex lastindex] = interpolate_IP_IA_and_IF(phi_unwrapped,r,t,Ts);
else
    firstindex = 1; 
    lastindex = size(ArcPoints,2);
    phi_unwrapped = zeros(1,length(phi));
    pause
end
ArcPoints = ArcPoints(:,firstindex:lastindex);
TangentPoints = TangentPoints(:,firstindex:lastindex);
phi_unwrapped = phi_unwrapped(firstindex:lastindex);
