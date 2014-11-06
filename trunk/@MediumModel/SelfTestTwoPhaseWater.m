function SelfTestTwoPhaseWater()
close all

%% compute liquid-vapour equlibrium using medium model
% note that some dilutant is required, otherwise the equilibrium will be
% pur liquid water between 0 and 100degC
water=MediumModel({'H2O','H2ObLb','N2'});
water.setT([linspace(0,200,25), linspace(99,101,25)]'+273.15);
water.setZ([0.495 0.495 0.01]);
nu=[-1 +1 0];
water.setNu(nu);
water.solveEq

%% Compute saturation temperature using XSteam
pVap=water.aeq(:,1)*water.P0*1e-5;
% condition <0.98 is needed because after this point the vapour is no longer saturated, so it is not a vlid comparison
swtValid=isfinite(pVap) & pVap<0.98; 
tSat_XSteam=nan(size(pVap));

for i=find(swtValid)'
    tSat_XSteam(i)=XSteam('Tsat_p',pVap(i));
end


figure
plot(water.T-273.15,water.aeq(:,1)*water.P0*1e-5)
hold all
plot(tSat_XSteam,water.aeq(:,1)*water.P0*1e-5)
xlabel('Temperature (degC)')
ylabel('Saturation vapour pressure (bar)')
grid on
legend({'MediumModel','XSteam'})

figure(2);clf;
plot(water.h,water.T-273.15,'x')
hline(100,'r','100 degC')
xlabel('Molar enthalpy (J/mol)')
ylabel('Temperature (degC)')
grid on
%legend({'MediumModel','XSteam'})


if max(abs(tSat_XSteam(swtValid)+273.15-water.T(swtValid))./water.T(swtValid)) > 5e-3
     error('MediumModel:SelfTestTwoPhaseWater:WrongSatTemp','Saturation temperature does not match that predicted by XSteam')
end

%%
disp('MediumModel.SelfTestTwoPhaseWater -- Test Passed')



end