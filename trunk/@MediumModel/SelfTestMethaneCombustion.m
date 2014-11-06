function SelfTestMethaneCombustion()
close all

%% compute liquid-vapour equlibrium using medium model
% note that some dilutant is required, otherwise the equilibrium will be
% pur liquid fuelandair between 0 and 100degC
dmolReac=1; %molar flow of reactants
fuelandair=MediumModel({'H2O','H2ObLb','N2','O2','CH4','CO2'});
fuelandair.setT([25]+273.15);
fuelandair.setZ([0 0 0.8 0.2 0.05 0]/1.05);%lambda 2
nu=[-1  1  0  0 0 0;...
     2  0  0 -2 -1 1];
fuelandair.notCondensed(2)=0; 
fuelandair.setNu(nu);
dmolFuel=dmolReac*fuelandair.Z(5);
hReacPerMolFuel=fuelandair.h*dmolReac/(dmolFuel); %total enthalpy flow per mole of CH4
dmReac=dmolReac*fuelandair.mm;
fuelandair.solveEq

% Check all CH4 consumed

if abs(fuelandair.Zeq(5))>1e-3
     error('MediumModel:SelfTestMethaneComb','Methane not consumed despite oxygen excess')
end

%% check HHV computation
% Shift all water to liquid, as this is the assumption in HHV calcs
fuelandair.setZ([ 0 sum(fuelandair.Zeq(1:2)) fuelandair.Zeq(3:6)]);
dmolProd=dmReac/fuelandair.mm;
hProdPerMolFuel=fuelandair.h*dmolReac/(dmolFuel); %total enthalpy flow per mole of CH4
mmCH4=12+4;
hHHVMediumModel=(-hProdPerMolFuel+hReacPerMolFuel)*1000/mmCH4;% HHV, J/kg
hHHVLit=55.5e6; %Ref: Haywood thermodynamic tables, CUP, 3rd Ed

if abs((hHHVMediumModel-hHHVLit)/hHHVLit)>5e-3
     error('MediumModel:SelfTestMethaneComb','HHV of methane combustion not correctly calculated')
end


disp('MediumModel.SelfTestMethaneCombustion -- Test Passed')



end