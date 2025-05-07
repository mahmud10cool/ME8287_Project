function[c, ceq] = evaluateConstraints(x)
%% CONSTRAINTFUNCTION computes the constraints. The input is an array x 
% containing the free variables. The outputs are the evaluated values of
% the inequality constraint c and equality constraint ceq. MMATLAB
% reference:https://www.mathworks.com/help/gads/gamultiobj.html#bvf79ug-nonlcon

    %% Assign variables
    dm = x(1); % [mm]
    delta = x(2); % [mm]
    dsy = x(3); % [mm]
    dst = x(4); % [mm]
    wst = x(5); % [mm]
    ast = x(6); % [degrees]
    
    %% Compute the inequality and equality constraints

    rs = 180; % stator radius [mm]
    dsp = 4; % stator tooth tip dimension [mm]
    dso = 2; % stator tooth tip thickness [mm]
    Omega = 10000; % rotational speed [RPM]
    
    rsi = rs - dsy - dst - dsp; % stator inner bore radius [mm]
    rm = rsi - delta; % magnet tip radius [mm]
    c1 = Omega*2*pi/60*rm/1000 - 175; % tip speed constraint
    % disp('c1')
    % disp(c1)
    
    Kcu = 0.4;
    dw=0.324861*0.0254*exp(-0.115942*18); % wire diameter [m]
    Sc = 0.25*pi*dw^2*1e6;
    alpha_c = pi/6;
    zQ = floor((Kcu/(2*Sc))*(alpha_c*(rs-dsy)*dst - wst*dst - 0.5*alpha_c*dst^2)); % number of turns
    % disp(z_Q)
    % Sc = 0.82355; % cross sectional area [mm^2]
    % zQ = round(Kcu/Sc/2*(pi/12*dst^2 + (pi/6*ri - wst)*dst)); % number of turns
    % 

    I_hat = 5*Sc*sqrt(2); % peak current [A]
    z_C = 4;
    k_w1 = 1; % fundamanetal winding factor
    r = rsi - delta/2; % airgap mid radius [mm]
    c2 = (3*zQ*I_hat*z_C*k_w1)/(pi*sqrt(2)*r*1e-3) - 80*1000; % rms current...
    % loading constraint
    % disp('c2')
    % disp(c2)
    c3 = (wst/2) - rsi*sin(deg2rad(ast/2)); % tooth geometry constraint
    %c3 =-1; % you have already given us a code for that

    % disp('rsi')
    % disp(rsi)
    % disp('wst')
    % disp(wst)
    % minima_ast = 2*asin(wst/2/rsi);
    % disp('minima ast')
    % disp(minima_ast)
    % 
    % ast_actual = ast*pi/180;
    % disp('actual ast')
    % disp(ast_actual)
    % disp('c3')
    % disp(c3)

    c4 = -(rs - dsy - dst - dsp - delta - dm); % valid geometry check 
    % disp('c4')
    % disp(c4)

    c = [c1;c2;c3;c4];
    
    ceq = checkToothTip(rs/1000, dso/1000, dsp/1000, dst/1000, wst/1000, dsy/1000, ast*pi/180);
    % disp('ceq')
    % disp(ceq)

end