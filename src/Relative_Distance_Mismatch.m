% relative distance mismatch

function d = Relative_Distance_Mismatch(x_temp, Hx_temp, x_previous, Hx_previous, y, Hy, y_previous, Hy_previous)

    v1 = [y - y_previous; Hy - Hy_previous];
    
    v2 = [(x_temp-x_previous)' (Hx_temp-Hx_previous)'];        % (need to add the difference between the circle and the data) create vector for mapped point of interest
    
    d = ((v1(1)-v2(:,1)).^2 + (v1(2)-v2(:,2)).^2).^0.5;       % need to minimise this guy to make sure mapping is optimal, check to see if the vectors are right!
    
end