function SelfTestElectrodePotentials
%%

[vCell,facThermoEff]=ComputeCellMetrics(25);

if abs(vCell-1.2289)>0.001
    error('MediumModel:SelfTestElectrodePotentials:wrongVoltage','Error in computation of standard electrode potential at 25degC ')
end

if abs(facThermoEff-0.83)>0.01
    error('MediumModel:SelfTestElectrodePotentials:wrongEff','Error in computation of thermodynamic FC efficiency at 25degC ')
end

[vCell,facThermoEff]=ComputeCellMetrics(600);


if abs(vCell-1.04)>0.01
    error('MediumModel:SelfTestElectrodePotentials:wrongVoltage','Error in computation of standard electrode potential at 600degC ')
end

if abs(facThermoEff-0.70)>0.01
    error('MediumModel:SelfTestElectrodePotentials:wrongEff','Error in computation of thermodynamic FC efficiency at 600degC ')
end

disp('MediumModel.SelfTestElectrodPotentials -- Test passed')


function [vCell,facThermoEff]=ComputeCellMetrics(tStk)

% Physical constants
facFaraday =   96.4853;
hSubstHhvH2 = -141.9103; % standard enthalpy change of combustion for H2 (J/mmol)
mmH2 =    2.0140;
%% set up feedstock

O2=MediumModel({'O2'});
O2.setZ([1]);
O2.setT(tStk+273.15);

% select arbitrary molar flow rate of hydrogen
dmolH2=1;

%%
hydrogen=MediumModel({'H2'});
hydrogen.setZ(1);
hydrogen.setT(tStk+273.15)

%% set up stack AOG  model
% 
if tStk>100
    aog=MediumModel({'H2O'});
else
    aog=MediumModel({'H2ObLb'});
end
aog.setZ([1]);
aog.setT(tStk+273.15);
dmolAog=dmolH2;

%%
dmolO2=0.5*dmolAog;
%% Compute the maximum useful work between stack inlet and stack outlet

QStackGibbsDel=(aog.mu*dmolAog - hydrogen.mu*dmolH2 - O2.mu*dmolO2);
qStackGibbsDelPerH2=QStackGibbsDel/dmolH2;
vCell=-qStackGibbsDelPerH2/(2*facFaraday*1000);

QStackHDel=(aog.h*dmolAog - hydrogen.h*dmolH2 - O2.h*dmolO2);
qStackHDelPerH2=QStackHDel/dmolH2;

facThermoEff=QStackGibbsDel/(hSubstHhvH2*1e3*dmolH2*mmH2);
