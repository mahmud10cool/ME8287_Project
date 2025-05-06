function import_Recoma35E(matName)
%IMPORT_Recoma33E Summary of this function goes here
%   Detailed explanation goes here

mu_r = 1.08;
Hc = 880000; %Coercivity of the magnet, units of A/m
Conductivity = 0.111; %Conductivity, units of MS/m
lam_fill = 1; %Lamination fill =1 (Magnets are not laminated)
nstr = 1; %No. of strands (FEMM Scripting Manual says this should be set to 1 for Magnets Ref: Pg 10)
mi_addmaterial(matName, mu_r, mu_r, Hc, 0, Conductivity, 0, 0, lam_fill,...
    0, 0, 0, nstr);

end
