function SelfTestBoudouard()
close all

%% Bouduard Carbon Formation Example
tExcel  = [ ...
    25,        127,        227,        327,        427,        527,        614, ...
    627,        727,        827,        927,       1027,       1127,       1227, ...
    ];

bouduard=MediumModel({'CO','CO2','Cbgrb'});
bouduard.setT(tExcel+273.15);
bouduard.setZ([1 0 0]);
nu=[2 -1 -1];
bouduard.setNu(nu);

GrMM=bouduard.gibbs;
bouduard.solveEq
bouduard.gibbsPlot
figure(1),
clf
plot(bouduard.T-273.15,bouduard.aeq,'DisplayName',{'CO','CO2','C'})
hold all

Xcoco2 = [ ...
    1          0
    1          0
    1          0
    0.998      0.002
    0.98       0.02
    0.89       0.11
    0.679      0.321
    0.637      0.363
    0.282      0.718
    0.078      0.922
    0.02       0.98
    0.006      0.994
    0.002      0.998
    0.001      0.999 ];
Gr  = [ ...
    112.8,       95.7,       78.9,       62.1,       45.3,       28.5,       13.9, ...
    11.8,         -5,      -21.8,      -38.6,      -55.4,      -72.2,        -89, ...
    ];


plot(tExcel,Xcoco2,':','DisplayName',{'CO2(excel)','CO(excel)'})
legend('Location','Best')


figure(2),clf
hold all
plot(tExcel,Gr,'o-','DisplayName','Excel')
plot(bouduard.T-273.15,GrMM/1e3,'o-','DisplayName','NASA')
xlabel('Temperature [degC]')
ylabel('Gibb''s Energy [kj/mol]')
legend('Location','Best')

%% atomic balance checks
dmolCOStart=1; %mol/s
% now use mass balance to compute final molar flow
dmReactants=dmolCOStart*28.0101; %g/s
dmolProducts=dmReactants./bouduard.mm;

rDmolOErr=(dmolCOStart-dmolProducts.*bouduard.Zeq(:,1)-2*dmolProducts.*bouduard.Zeq(:,2))/dmolCOStart;

if rDmolOErr>1e-9
    error('MediumModel:SelfTestBoudouard:MolarBalance','Oxygen molar balance error')
end
    

%% Gibbs energy and composition checks


if max(abs(GrMM'/1e3 - Gr)./Gr)>0.1 % accept about 5% error in calc of gibbs energy
    error('MediumModel:SelfTestBoudouard:WrongGibbs','Gibbs free energy change of reaction does not match literature values')
end

if max(abs(Xcoco2-bouduard.aeq(:,[2 1])))>0.05 % accept about 5% error in calc of composition 
    error('MediumModel:SelfTestBoudouard:WrongComposition','Composition does not match literature values')
end
disp('MediumModel.SelfTestBoudouard -- Test Passed')

end