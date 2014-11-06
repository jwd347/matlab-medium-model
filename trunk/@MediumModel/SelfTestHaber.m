function SelfTestHaber()
close all

%% Model the Haber Process to compare with published kc
haber=MediumModel({'N2','H2','NH3'});
haber.setT(([300 400:50:600]')+273.15);
nu=[-1 -3 2];
haber.setNu(nu);

haber.gibbs;
haber.solveEq;
kc = exp(haber.ln_kc) ;

% Published in
% Chemistry the Central Science" Ninth Ed., by: Brown, Lemay, Bursten, 2003, ISBN 0-13-038168-3
kcVal = [4.34e-3 1.64e-4 4.51e-5 1.45e-5 5.38e-6 2.25e-6];
figure
semilogy(haber.T-273.15,kcVal,'o-',haber.T-273.15,kc,':.')
legend('Published Value','Predicted Value')
ylabel('K_c')
xlabel('Temperature [degC]')
title('Comparison of Equilibrium constants for the Haber process')

if max(abs(haber.ln_kc - log(kcVal')))>0.15 % accept about 16% error in calc of kc
    error('MediumModel:SelfTestHarber:WrongLnKc','kc for Harber process does not match values from literature')
end

%% Check ATE calc functionality
tCatalyst=500;
[tEq,tATE]=equilibriumTemperature(haber,haber.Zeq,tCatalyst);

if max(abs(tEq-haber.T))>1
    error('MediumModel:SelfTestHarber:ATECalc','Failure to correctly compute equilibrium temperature from a given composition')
end

if max(abs(tEq+tATE-tCatalyst))>1e-6
    error('MediumModel:SelfTestHarber:ATECalc','Failure to correctly compute ATE from a given composition')
end


disp('MediumModel.SelfTestHaber -- Test Passed')