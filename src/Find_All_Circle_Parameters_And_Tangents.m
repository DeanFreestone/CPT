
function [x0_temp Hx0_temp r_temp phi_temp x_temp Hx_temp M_temp e d Norm] ...
    =  Find_All_Circle_Parameters_And_Tangents(MaxPoints, InitialPoints, PointsStep, ...
    y_s_all, Hy_s_all, x_previous, Hx_previous)

    center = floor(length(y_s_all)/2) + 1;              % length of y_s_all is odd so there is a centre

    PointsFromCenterOuter = floor(MaxPoints/2);         % should be m_{max}
    PointsFromCenterInner = floor(InitialPoints/2);         % m_{min}

    NPoints = floor((PointsFromCenterOuter - PointsFromCenterInner)/PointsStep) + 1;          % how many combos we have to try

    DistanceFromCenter = PointsFromCenterInner:PointsStep:PointsFromCenterOuter;
        
    x0_temp = zeros(1,NPoints);                   % initialise to zeros for speed. 
    Hx0_temp = zeros(1,NPoints);
    r_temp = zeros(1,NPoints);
    phi_temp = zeros(1,NPoints);
    x_temp = zeros(1,NPoints);
    Hx_temp = zeros(1,NPoints);
    e = zeros(1,NPoints);
    d = zeros(1,NPoints);
    M_temp = zeros(1,NPoints);
    Norm = zeros(2,NPoints);
    
    PointSpace = 1;
    
    for n=1:NPoints         %PointsFromCenterInner:PointsStep:PointsFromCenterOuter
        
        StartSegment = center-DistanceFromCenter(n);
        EndSegment = center+DistanceFromCenter(n);
        
        y_s = y_s_all(StartSegment:PointSpace:EndSegment);                 % put the y and Hy signal in a more readable form
        Hy_s = Hy_s_all(StartSegment:PointSpace:EndSegment);               % this is the segment used to fit the circle

        Tangent = [y_s(end) - y_s(1) Hy_s_all(end) - Hy_s(1)];
        Norm(:,n) = [Tangent(2); -Tangent(1)];
        
        M_temp(n) = 2*DistanceFromCenter(n)+1;          % size of segment
        
        XY = [y_s Hy_s];                                                % segment in matrix so we can use the circle fit

        Par = CircleFitByTaubin(XY);                    % Nonlinear. from Mathworks online. Uses LMS fit

        x0_temp(n) = Par(1);                                % estimate for y0 translation
        Hx0_temp(n) = Par(2);                               % estimate for Hy0 translation
        r_temp(n) = Par(3);                                 % estimate for radius (envelope)
        
        % might be able to move these guys outside the loop to make it
        % faster (to do later)
%         ~~~~~~~~~~~
        phi_temp(n) = atan2(Hy_s_all(center)-Hx0_temp(n), y_s_all(center)-x0_temp(n));
        x_temp(n) = r_temp(n)*cos(phi_temp(n));
        Hx_temp(n) = r_temp(n)*sin(phi_temp(n));
        
        % maybe we should use the edge of the circle for this guy
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        e(n) = mean( (y_s_all([StartSegment,center,EndSegment]) - x0_temp(n)).^2 ...
            + (y_s_all([StartSegment,center,EndSegment]) - Hx0_temp(n)).^2 - r_temp(n)^2 );   % estimate the mse for the signal
        
% new bit for 130032010
% ~~~~~~~~~~~~~
% for the RDM
% ~~~~~~~~~~~~~~~~
      
        x_previous_s = x_previous(end-DistanceFromCenter(n)+1:PointSpace:end);           % take the part of the previous estimate that correspond to the segment size
        Hx_previous_s = Hx_previous(end-DistanceFromCenter(n)+1:PointSpace:end);
        
        % here we remove all the previous values that were nans
        x_previous_s_not_nan = x_previous_s(~isnan(x_previous_s));          % cant use any values that are NaN 
        
        if ~isempty(x_previous_s_not_nan)
            
            Hx_previous_s_not_nan = Hx_previous_s(~isnan(Hx_previous_s));   % need to use the value as close the the first sample that is not NaN
        
            y_previous_s = y_s_all(center-DistanceFromCenter(n):PointSpace:center-1);      % cant include current sample!
            Hy_previous_s = Hy_s_all(center-DistanceFromCenter(n):PointSpace:center-1);

            y_previous_s_not_nan = y_previous_s(~isnan(x_previous_s));          % make sure we are dealing with the same samples when we make our comparisons for the RDM 
            Hy_previous_s_not_nan = Hy_previous_s(~isnan(Hx_previous_s));  
            
            % changing this from d_temp to d
            d(n) = Relative_Distance_Mismatch(x_temp(n), Hx_temp(n), x_previous_s_not_nan(1), ...
                Hx_previous_s_not_nan(1), y_s_all(center), Hy_s_all(center), y_previous_s_not_nan(1), ...
                Hy_previous_s_not_nan(1));

%             for nn=1:length(x_previous_s_not_nan)
%                 d_temp(nn) = Relative_Distance_Mismatch(x_temp(n), Hx_temp(n), x_previous_s_not_nan(nn), Hx_previous_s_not_nan(nn), y_s_all(center), Hy_s_all(center), y_previous_s_not_nan(nn), Hy_previous_s_not_nan(nn));
% %             d(m) = Relative_Distance_Mismatch(x_temp(m), Hx_temp(m), x_previous_s_not_nan(end), Hx_previous_s_not_nan(end), y_s_all(center), Hy_s_all(center), y_previous_s_not_nan(end), Hy_previous_s_not_nan(end));
%             end
            
%             d(n) = mean(d_temp);
        else
            d(n) = nan;
        end
%         ~~~~~~~~~~
%         m = m+1;
    end
    
%     d=e;
end