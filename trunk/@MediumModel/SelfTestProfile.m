function SelfTestProfile()
close all


%% Reformate Composition
%%%  Define the mixture model for reformate
% Use this vector order always H2, CH4, C0, CO2, H2O
fuel=MediumModel({'H2','CH4','CO','CO2','H2O','N2'});
% Assume Steam:Carbon=2.5
Z=[0 1 0 0 2.5 0.002];
fuel.setZ(Z./sum(Z));
fuel.setT([500:1:700]'+273.15);

nu=[    [3 -1 1  0 -1 0 ]; ...
        [1 0  -1 1 -1 0]];
    
fuel.setNu(nu);
fuel.gibbs;

tiStart=tic;
fuel.solveEq
tiSolveEq=toc(tiStart)

tiStart=tic;
fuel.solveEqQSearch;
tiSolveEqQSearch=toc(tiStart)



disp('MediumModel.SelfTestProfile -- Test Executed')
