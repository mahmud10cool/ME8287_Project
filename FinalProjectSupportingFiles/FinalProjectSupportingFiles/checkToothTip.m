function [ceq] = checkToothTip(statorRadius, dso, dsp, dst, wst, dsy, ast)
%checkToothTip Validates the tooth tip geometry
%   Returns 0 if the tooth tip geometry is valid (thinnest portion of tooth
%   tip is thicker than dso/2
%   Returns 1 if the tooth tip geometry is not valid.
%   Units:
%       - statorRadius, dso, dst, wst, dsy are all in [m].
%       - ast is the angle alpha_st in [radians].

rs = statorRadius; 
rsi = rs - dsy - dst - dsp;

ceq = 1;

%only run the test if wst is valid.
if wst < ast*rsi 
    x0 = wst/2;
    astp = asin(x0/(rsi+dsp));
    y0 = (rsi+dsp)*cos(astp);

    x2 = rsi*sin(ast/2);
    y2 = rsi*cos(ast/2);

    x1 = x2;
    y1 = y2 + dso;

    m = (y1 - y0)/(x1 - x0);
    b = y1 - m*x1;

    x = -m*b/(1+m^2);    

    rtx = sqrt( (1+m^2)*x^2 + b^2 + 2*m*b*x);    %radius at minimum x coord
    rtx1 = sqrt( (1+m^2)*x1^2 + b^2 + 2*m*b*x1); %radius at tooth edge
    if ((x > x1) && (rtx1 - rsi > dso/2)) || ((x < x1) && (rtx - rsi > dso/2))
        %We passed the test!            
        ceq = 0;    
    else
        %We failed the test!
        ceq = 1;
    end
    
    %debug info: 
    %fprintf('---> ceq=%g, ast=%g, astp=%g, x = %g, x1=%g, x0=%g, y1=%g, y0=%g, r(x)=%g, r(x1) - rsi = %g, r(x) - rsi = %g \n', ceq, ast*180/pi, astp*180/pi, x, x1, x0, y1, y0, rtx, rtx1-rsi, rtx-rsi)
end
end