function SelfTestViscTCondRho
close all
% Tests of Dynamic viscosity, Thermal conductivity & Gas density

%% Test point data
% from google
rhoAir_iupac = 1.2754; % @ 1bar, 0C
rhoAir_nist = 1.2922; % @ 1atm, 0C
rhoAir_isa = 1.225; % @ 1atm, 15C
% from other
rhoAir_300C = 0.6158;
rhoAir_600C = 0.4042;
rhoAir_80kPa = 0.95;

% from http://www.lmnoeng.com/Flow/GasViscosity.php
dViscAir_0C = 0.0173*1e-3; 
dViscAir_300C = 0.0299*1e-3; 
dViscAir_600C = 0.0393*1e-3;

% from http://www.nist.gov/data/PDFfiles/jpcrd269.pdf
kAir_300K = 26.19*1e-3; 
kAir_600K = 45.69*1e-3; 
kAir_1000K = 67.85*1e-3; 
 
%% MediumModel tests

% rho of gases
air = MediumModel({'N2','O2','Ar','CO2'});
air.setTandZ(273.15, [0.78084 0.20944 0.00933 0.00039]); % Top 4 constituents
assert(abs(air.rho-rhoAir_iupac)<1e-5,'Wrong Air density at IUPAC condition.');

air.setP(101325);
assert(abs(air.rho-rhoAir_nist)<1e-3,'Wrong Air density at NIST condition.');

air.setT(288.15);
assert(abs(air.rho-rhoAir_isa)<1e-4,'Wrong Air density at ISA condition.');

air.setP(101325);
air.setT(300+273.15);
assert(abs(air.rho-rhoAir_300C)<1e-4,'Wrong Air density at 300C condition.');

air.setT(600+273.15);
assert(abs(air.rho-rhoAir_600C)<1e-4,'Wrong Air density at 600C condition.');

air.setP(80000);
air.setT(20+273.15);
assert(abs(air.rho-rhoAir_80kPa)<1e-3,'Wrong Air density at 80kPa,20C condition.');

% dynVisc of gases
air.setP(100000);
air.setT(273.15);
assert(abs(air.dynVisc-dViscAir_0C)<1e-6,'Wrong Air viscosity at 0C condition.');
air.setT(300+273.15);
assert(abs(air.dynVisc-dViscAir_300C)<1e-6,'Wrong Air viscosity at 300C condition.');
air.setT(600+273.15);
assert(abs(air.dynVisc-dViscAir_600C)<1e-6,'Wrong Air viscosity at 600C condition.');

% k of gases
air.setP(100000);
air.setT(300);
assert(abs(air.k-kAir_300K)<1e-3,'Wrong Air thermal conductivity at 300K condition.');

air.setT(600);
assert(abs(air.k-kAir_600K)<1e-3,'Wrong Air thermal conductivity at 600K condition.');

air.setT(1000);
assert(abs(air.k-kAir_1000K)<1e-3,'Wrong Air thermal conductivity at 1000K condition.');

%% rho,dynVisc,k of gases & liquid
air_water = air.copy;
air_water.addSpecies({'H2ObLb'});
Z_dryair = [0.78084 0.20944 0.00933 0.00039];
air_water.setZ([Z_dryair 0.5]/(sum(Z_dryair)+0.5));
assert(abs(air.rho-air_water.rho)<1e-10,'Wrong rho calc');
assert(abs(air.dynVisc-air_water.dynVisc)<1e-10,'Wrong dynVisc calc');
assert(abs(air.k-air_water.k)<1e-10,'Wrong k calc');

% Tests of mixtures with gases having no dynVisc,k
gas = MediumModel({'N2','NO'});
gas.setTandZ(300,[0.8 0.2]);
% dynVisc of gases
assert(isnan(gas.dynVisc),'gas dynVisc expected as NaN')
% k of gases
assert(isnan(gas.k),'gas k expected as NaN')
% rho of gases
assert(abs(gas.rho-1.1391)<1e-4,'gas rho expected as 1.1391')

% Tests with vectors
air.setT((300:100:600)');
air.rho
air.dynVisc
air.k

air.setTandZ((300:100:600)',[0.7 0.2 0.05 0.05;...
                             0.6 0.2 0.1 0.1;...
                             0.5 0.2 0.15 0.15;...
                             0.6 0.1 0.2 0.1;]);
air.rho
air.dynVisc
air.k

disp('MediumModel.SelfTestViscTCondRho -- Test completed')
