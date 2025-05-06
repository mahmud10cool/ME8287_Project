function [designEval] = evaluateDesign(materials, dimensions, p, winding, settings)

% Co-author: Mahmud Suhaimi Ibrahim
% Date: 4/10/2025

%EVALUATEDESIGN Evaluates a candidate motor design.
% materials is a structure containing the material names for the magnet, coil and iron.
% dimensions is a structure containing the fimensions of the machine in
% units of [mm]. 
% p is the number of pole pairs
% winding is a structure that contains information about the stator
% winding configuration. 
% settings is a structure that contains values pertaining to the current 
% excitation applied to the stator coils and the motor speed in RPM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the groups to which each of the geometry belongs. This will help
% compute cross-section area easily using FEMM post processing function
rotorIronGroup = 10; %Group containing the rotor iron == 10
statorIronGroup = 1; %Group containing the statorIron == 1
slotGroup = 2; % Group containing the winding labels == 2 

% helpful unit definitions
PI=pi; Pi=pi;
deg=Pi/180; % degree to radian

% Find the No. of slots (Q) Hint: use the winding structure to do
% this
Q = length(winding.topSlots.circuit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Initialize FEMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open FEMM, Create a new magnetics problem, define length units and type
% of problem (planar) 
openfemm()
newdocument(0);
mi_probdef(0,'millimeters','planar',1e-8,dimensions.length,-30,0);
smartmesh(0); % Turns off the smart mesh

%% Fetch Library materials
mi_getmaterial('Air');
mi_getmaterial(materials.coil.name);
import_Recoma35E(materials.magnet.name);
import_M19_29Ga(materials.iron.name);   %Imports the BH curve and adds the material

%% Create Circuits
mi_addcircprop('U', 0, 1); %Let initial currents be 0. We will assign them later
mi_addcircprop('V', 0, 1);
mi_addcircprop('W', 0, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Step (a): Create Stator %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makeStator(materials, dimensions.stator, dimensions.outerRadius,Q, winding);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Step (b): Create Rotor  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimensions.rotor.OuterRadius = dimensions.outerRadius - (dimensions.stator.dsp + dimensions.stator.dst + dimensions.stator.dsy + dimensions.delta);
makePMRotor(materials, dimensions.rotor,dimensions.rotor.OuterRadius, p, rotorIronGroup)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Step (c): Create Air regions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add air to the inner bore of the rotor
mi_addblocklabel(0,0);
mi_selectlabel(0,0);
mi_setblockprop('Air', 1, 0, 0, 0, 1, 0);
mi_clearselected;

% Add air to the airgap
mi_addblocklabel(dimensions.rotor.OuterRadius+0.5,0);
mi_selectlabel(dimensions.rotor.OuterRadius+0.5,0);
mi_setblockprop('Air', 1, 0, 0, 0, 1, 0);
mi_clearselected;

% Add air just outside stator
mi_addblocklabel(dimensions.outerRadius+0.5,0);
mi_selectlabel(dimensions.outerRadius+0.5,0);
mi_setblockprop('Air', 1, 0, 0, 0, 1, 0);
mi_clearselected;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Step (d): Set Boundary Condition %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi_makeABC(7,dimensions.outerRadius*1.5,0,0,1);
mi_saveas('femmTest1.fem');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Step (e): Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the initial parameters for the FEA solves
n = 360/settings.lowestHarmonic;  % angle through which to spin the rotor in degrees
steps = settings.steps; % number of steps the rotor is evaluated at
deltaThetaMec = (360/settings.lowestHarmonic)/steps; % angle increment of theta_mec in degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform a series of finite element analyses
% First Run an analysis through an entire spin of the rotor and record
% element centroid flux density and vector potential (using the mesh from the first iteration)
% at every step.  This information will then be used to estimate losses
opendocument('femmTest1.fem');

mi_saveas('temp.fem');

for thetaMec = 0:deltaThetaMec:(n-1) % vary the rotor position
    kk = round(thetaMec/deltaThetaMec + 1);
    
    setCurrents(settings.iHat, p*thetaMec, settings.phi);  % Set currents
    
    % Mesh and solve the model (add code)
    mi_analyze(1);   
    % Load solution (add code)
    mi_loadsolution;

    if (thetaMec == 0)
        nn = mo_numelements;  % Record the initial mesh elements if the first time through the loop
        B = zeros(floor(n/deltaThetaMec),nn); % define a matrix that will hold the flux density info
        A = zeros(floor(n/deltaThetaMec),nn); % define a matrix that will hold the vector potential info
        centroid = zeros(nn,1); % define a matrix that will hold the centroid coordinates
        a = zeros(nn,1); % define a matrix that will hold the mesh element area
        group_num = zeros(nn,1); % define a matrix that will hold the group number the mesh element belongs to
        torque_array = zeros(floor(n/deltaThetaMec),1);

        for m = 1:nn
            elm = mo_getelement(m);
            % centroid is a vector of complex numbers that represents the location of
            % the centroid of each element.  The real part is x, the
            % imaginary part is y.  The purpose of representing the
            % location in this way is that it is easy to rotate the
            % location of the centroid (by multiplying by a complex number)
            % to find the point that corresponds to the centroid when the
            % rotor is rotated.
            centroid(m) = elm(4) + 1j*elm(5);
            a(m) = elm(6); % element area in the units used to draw the geometry
            group_num(m) = elm(7); % group number associated with the element
        end
    end

   % Store element flux densities
   u=exp(1j*thetaMec*deg);
   for m = 1:nn
       % Element is on the rotor if group number is >=rotorIronGroup.   
       if(group_num(m)>=rotorIronGroup)
            % for rotor elements, the location in the original mesh is rotated so that the
            % flux density at the same point with the rotor rotated through
            % angle k can be evaluated. Multiplying by the complex number u
            % does this rotation. 
            p_2 = centroid(m)*u;
    
            % Flux densities bx and by are evaluated and rolled into a complex number.
            % Dividing by the complex number u rotates the flux density
	        pv = mo_getpointvalues(real(p_2),imag(p_2));
    
            if (group_num(m)==rotorIronGroup)
                % store flux density for elements in rotor core
		        B(kk,m) = (pv(2)+1j*pv(3))/u;
            else
                % store vector potential for elements that are in PMs
		        A(kk,m) = pv(1);
            end
        elseif (group_num(m) > 0) 
            % element is on the stator
            % since the stator doesn't move, there is no need for the
            % reference frame transformations
            p_2 = centroid(m);
            B(kk,m) = (mo_getb(real(p_2),imag(p_2))*[1;1j]);
       end       
    end 

    %get the problem info
    probinfo=mo_getprobleminfo;

	% select all blocks on the rotor and compute the torque
	for block=rotorIronGroup:(rotorIronGroup+2*p)
        mo_groupselectblock(block);
    end
    torque_array(kk)=mo_blockintegral(22);  
    mo_close;

    % rotate the rotor to the position for the next iteration.
    for block=rotorIronGroup:(rotorIronGroup+2*p)
        mi_selectgroup(block);
    end
    mi_moverotate(0, 0, deltaThetaMec); 
	mi_clearselected();
end

closefemm;   % close FEMM after finite element runs are finished
delete('temp.fem');
delete('temp.ans');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Store Data for post processing %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feaData.A = A; % Vector potential
feaData.B = B; % Flux density
feaData.L = probinfo(3); %Axial length [m]
feaData.units = probinfo(4); % length units
feaData.eleArea = a; %element area in [m^2]
feaData.groupNum = group_num; % group numbers
feaData.centroid = centroid; % centroid of mesh elements
save('FEA_Data.mat'); % Save the FEA Data. This helps you test and debug the 
                      % post-processing part of the code without requiring
                      % to re-run FEA again - saves a lot of time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Post Processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Compute avg. torque and the torque ripple%%%%%%%%%%%%%%%%%%%%%
designEval.torque.average = mean(torque_array);
designEval.torque.ripple = (abs(max(torque_array)-min(torque_array))/designEval.torque.average);

%%%%%%%%%%%%%%%%%%%%%%% Compute Core Losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the square of the amplitude of each harmonic at the centroid of
% each element in the mesh. Matlab's built-in FFT function makes this easy.
% refer the FEMM example for this
ns=n/deltaThetaMec;

Bxfft=abs(fft(real(B)))*(2/ns);
Byfft=abs(fft(imag(B)))*(2/ns);
Bsq=(Bxfft.*Bxfft) + (Byfft.*Byfft);

% Compute the volume of each element in units of meter^3
h = probinfo(3);            % Length of the machine in the into-the-page direction
lengthunits = probinfo(4);  % Length of drawing unit in meters
v = a*h*lengthunits^2;

%%%%%%%%%%%%%%%%%%%%%%% Compute Magnet Losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute fft of A at the center of each element
Jm=fft(A)*(2/ns);
for k=1:2*p
	g3=(group_num==(10+k));
	% total volume of the magnet under consideration;
	vmag=v'*g3;
	% average current in the magnet for each harmonic
	Jo=(Jm*(v.*g3))/vmag;    
	% subtract averages off of each each element in the magnet
	Jm = Jm - Jo*g3';
end

%%%%%%%%%%%%%%%%%%%%%%%%% Compute copper loss%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dwire=0.324861*0.0254*exp(-0.115942*materials.coil.AWG); % wire diameter in meters as a function of AWG
Acond = 0.25*pi*dwire^2;
owire = materials.coil.conductivity; % conductivity of the wire in S/m at prescribed deltaT
zQ = abs(winding.topSlots.zQ(1));
kov = 1.8;
r_slot = dimensions.outerRadius - dimensions.stator.dsy - 0.5*dimensions.stator.dst;
alpha_c = (2*pi)/Q;
Tau_u = alpha_c * r_slot; % in [mm]
l_c = 1e-3*(2*dimensions.length + pi*((dimensions.stator.wst + Tau_u)/2) + 2*kov*Tau_u*(winding.y-1)); % in [m]

r_coil = (zQ*l_c)/(owire*Acond);
Iphase=settings.iHat/sqrt(2);
winding_loss = Q*(r_coil)*Iphase^2;

ch = materials.iron.ch;
ce = materials.iron.ce;
sf = materials.iron.sf;

omag = materials.magnet.conductivity;

rho = materials.iron.density;

% results=[];

items = 37;

range_of_speed = linspace(0,settings.RPM,items);

total_loss = NaN(size(range_of_speed));
outputPower = NaN(size(range_of_speed));
efficiency = NaN(size(range_of_speed));
rotor_loss = NaN(size(range_of_speed));
magnet_loss = NaN(size(range_of_speed));
stator_loss = NaN(size(range_of_speed));

for i = 1:length(range_of_speed)

    thisSpeed = range_of_speed(i);
    thisFrequency = thisSpeed/60; % mechanical speed in Hz
    
    w=0:(ns-1);
    w=settings.lowestHarmonic*thisFrequency*w.*(w<(ns/2));  
    
    % Now, total core loss can be computed in one fell swoop...
    % Dividing the result by sf corrects for the lamination stacking factor
    g1=(group_num==10);
    rotor_loss(i) = ((ch*w+ce*w.*w)*Bsq*(rho*v.*g1))/sf;
    
    g2=(group_num==1);
    stator_loss(i) = ((ch*w+ce*w.*w)*Bsq*(rho*v.*g2))/sf;
    
    % Add up eddy current losses in the magnets
    magnet_loss(i) = (1/2)*((omag*(2*pi*w).^2)*(abs(Jm).^2)*v);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Compute total loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    total_loss(i) = rotor_loss(i) + stator_loss(i) + winding_loss + magnet_loss(i);
    
    outputPower(i) = designEval.torque.average*thisSpeed*(2*pi/60);
    efficiency(i) = 100*(outputPower(i)/(outputPower(i)+total_loss(i)));

    % results = [results; thisSpeed, rotor_loss, stator_loss, magnet_loss, winding_loss, total_loss];
end

save efficiency_data.mat total_loss range_of_speed efficiency outputPower rotor_loss magnet_loss stator_loss winding_loss

designEval.loss.winding = winding_loss;
designEval.loss.magnets = magnet_loss;
designEval.loss.rotorIron = rotor_loss;
designEval.loss.statorIron = stator_loss;

%%%%%%%%%%%%%%%%%%%%%Efficiency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
designEval.efficiency = efficiency(end);

end