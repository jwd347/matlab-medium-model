%% Medium Model Tutorial
%%
% Work in progress. Contains LaTeX tags for publishing.
function Tutorial()

%% Reformer model
% \subsection{Parameters}


%% 
% \subsection{Reformate}
% General approach is to setup mixture with correct S:C and then equilibrate at tFuelRefOut 
% 
% First define some constants and some flow rates :

PhysConst; % load physical constants
facStToC =2.8; % S:C set-point
tFuelRefOut=[500 550 600 650]; % vector of temperatures at which to equilibrate the mixture
dmolFuelRef=1;%arbitrary flow rate. 
dmolWatFc= dmolFuelRef*facStToC;
%% 
% Create an instance of the medium model containing the chemical species of
% interest:
reformate=MediumModel({'H2','CH4','CO','CO2','H2O'});
%%
% Define the stochiometry matrix for WGS and reformation reactions:
nu=[    [3 -1 1  0 -1  ]' ...
    [1 0  -1 1 -1 ]'];
% This is an NxM matrix where M is the number of species and N is the
% number of reactions. The entry in row $n$ and column $m$ indicates the
% number of molecules of the n$^{th}$ species produced by the $m^{th}$
% reaction.
reformate.setNu(nu);
reformate.setT(tFuelRefOut+273.15);
reformate.setZ([0 1 0 0 facStToC ]'./(1+facStToC));
reformate.solveEq;
%%
% Having 'solved' the medium model, there are a number of properties than
% may be examined:
reformate.h % molar enthalpy, J/mol
reformate.mu % Chemical potential, J/mol
reformate.Zeq % equilibrium composition
%%
% A full list of properties are available by typing "reformate." and then
% TAB at the command prompt. We can use the above calcs to compute reformer
% yield as a function of temperature:
figure(1);clf
%plot(....
%%
% Perform a mass balance to compute the molar flow rate:
dmReformate= (dmolFuelRef*mmCH4+dmolWatFc*mmH2O)*1e-3;
dmolReformate=dmReformate*1e3./reformate.mm;





%% Methane combustion




