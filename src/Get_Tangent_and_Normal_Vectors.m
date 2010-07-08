% Samples_For_Tangent should be an odd number such that it has a center

% Note, When using arcs that are semi-circles:
% we could simply use a line between samples end samples of the arc but
% this does not allow for noise or disturbance.
% a simple regression line with vertical offset does not work when the
% tangent is vertical.
% a regression line with perpendicular offsets when the tangent is horizontal 
% and vertical, but not for 45 degree tangents, and we have to choose
% between two solulations.

% seems the best solutions is to simply use two points in the segment and
% make a vector that runs through them. If they are far enough apart noise
% should not be too much off a problem.

function [Simple_Tangent Simple_Norm] = Get_Tangent_and_Normal_Vectors(Samples_For_Tangent,y_s_all, Hy_s_all)

    center = floor(length(y_s_all)/2) + 1;              % length of y_s_all should be odd so there is a centre
    
    Samples_Either_Side = floor(Samples_For_Tangent/2);
    
%     y_t = y_s_all(center - Samples_Either_Side:center + Samples_Either_Side);
%     Hy_t = Hy_s_all(center - Samples_Either_Side:center + Samples_Either_Side);
% 
%     b_ls = regress(Hy_t',[ones(size(y_t))' y_t']);            % y=ax+b or y = p(1) + p(2)x
%     tangent = b_ls(2)*y_t + b_ls(1);
% 
%     % find the normal vector (z) of the estimated tangent
%     z = [y_t(end)-y_t(1) tangent(end)-tangent(1)];
%     a = [z(2); -z(1)];                             % (normal to the separating line) translated to the origin
    
    
    Simple_Tangent = [y_s_all(center + Samples_Either_Side) - y_s_all(center - Samples_Either_Side)...
    Hy_s_all(center + Samples_Either_Side) - Hy_s_all(center - Samples_Either_Side)];
    
    Simple_Norm = [Simple_Tangent(2); -Simple_Tangent(1)];

%     Simple_Tangent = Simple_Tangent/norm(Simple_Tangent);
%     Simple_Norm = Simple_Norm/norm(Simple_Norm);
%     
%     z = Simple_Tangent;
%     a = Simple_Norm;

end
