
function lambda = Check_Circle_Center(y_n, Hy_n, x0_temp, Hx0_temp, x_temp, Hx_temp, a)

    Trans_x0 = x0_temp - y_n;       % the mapped offset translated or reference to origin (since normal vector is at origin).
    Trans_Hx0 = Hx0_temp - Hy_n;      

    % now we check the angle between the normal vector to the fitted
    % line and mapped circle centre  
    
%     for n=1:length(x0_temp)
%         lambda(n) = -a(:,n)'*[x_temp(n); Hx_temp(n)];
%     end
    
%     X = [Trans_x0' Trans_Hx0'];       % should be a mx2 matrix, with Trans_x0 filling a column
%     lambda_y = (X*a)';                           % a is a 2x1 vector, which is the normal vector to the tangent of the local arc
%     
    X = -[x_temp; Hx_temp]';
    lambda = (X*a)';
    
end