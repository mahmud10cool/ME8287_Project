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

    %% Constraint 1
    rs = 180; % stator radius [mm]
    dsp = 4; % stator tooth tip dimension [mm]
    dso = 2; % stator tooth tip thickness [mm]
    Omega = 10000; % rotational speed [RPM]
    
    rsi = rs - dsy - dst - dsp; % stator inner bore radius [mm]
    rm = rsi - delta; % magnet tip radius [mm]
    c1 = Omega*((2*pi)/60)*(rm*1e-3) - 175; % tip speed constraint

    %% Constraint 2
    Kcu = 0.4;
    dw=0.324861*0.0254*exp(-0.115942*18); % wire diameter [m]
    Sc = 0.25*pi*dw^2*1e6;
    alpha_c = pi/6;
    zQ = floor((Kcu/(2*Sc))*(alpha_c*(rs-dsy)*dst - wst*dst - 0.5*alpha_c*dst^2)); % number of turns


    I_hat = 5*Sc*sqrt(2); % peak current [A]
    z_C = 4;
    k_w1 = 1; % fundamanetal winding factor
    r = rsi - delta/2; % airgap mid radius [mm]

    % RMS current
    c2 = (3*zQ*I_hat*z_C*k_w1)/(pi*sqrt(2)*r*1e-3) - 80*1e3; 

    %% Constraint 3
    % Tooth geometry constraint
    c3 = (wst/2) - rsi*sin(deg2rad(ast/2));

    %% Constraint 4
    % Valid geometry check
    c4 = -(rs - dsy - dst - dsp - delta - dm); 

    c = [c1;c2;c3;c4];
    
    %% Constraint 5
    ceq = checkToothTip(rs/1000, dso/1000, dsp/1000, dst/1000, wst/1000, dsy/1000, ast*(pi/180));

end