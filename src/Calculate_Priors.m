

function [p1 p2] = Calculate_Priors(r_temp, phi_temp, r_previous, phi_previous)

    dr = r_temp - r_previous;
    p1 = (phi_temp+pi) - (phi_previous+pi);  % check to see if the phase is increasing (obviously will not work at point of 2pi phase transition)

    p2 = abs(p1) - 1*abs(dr);        % motivated by the Bedrosian theorem, must be > 0
    
end
    
 