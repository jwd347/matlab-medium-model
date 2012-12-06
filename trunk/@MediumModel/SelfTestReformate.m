function fuel=SelfTestReformate()
close all


%% Reformate Composition


%%% Reformer Data From HSC

%Temperature	H2(g)	CH4(g)	CO(g)	CO2(g)	N2(g)
%C	mol-%	mol-%	mol-%	mol-%	mol-%

%    H2(g)	H2O(g)	CH4(g)	CO(g)	CO2(g)	N2(g)

hscData = [
450	25.730	49.260	18.230	0.424	6.115	0.2476
455	26.590	48.530	17.870	0.478	6.289	0.2463
460	27.450	47.800	17.510	0.537	6.461	0.2450
465	28.330	47.060	17.140	0.602	6.630	0.2437
470	29.210	46.320	16.760	0.673	6.796	0.2423
475	30.090	45.570	16.390	0.752	6.959	0.2410
480	30.980	44.820	16.000	0.838	7.117	0.2396
485	31.870	44.070	15.620	0.932	7.270	0.2382
490	32.770	43.320	15.220	1.033	7.417	0.2367
495	33.670	42.560	14.830	1.144	7.559	0.2353
500	34.570	41.810	14.430	1.264	7.694	0.2339
505	35.470	41.060	14.020	1.394	7.822	0.2324
510	36.370	40.310	13.610	1.533	7.943	0.2309
515	37.270	39.560	13.200	1.683	8.055	0.2294
520	38.160	38.820	12.790	1.844	8.158	0.2279
525	39.060	38.080	12.370	2.015	8.253	0.2264
530	39.940	37.340	11.950	2.198	8.337	0.2249
535	40.830	36.620	11.530	2.392	8.412	0.2233
540	41.700	35.900	11.110	2.598	8.477	0.2218
545	42.570	35.190	10.680	2.814	8.531	0.2203
550	43.420	34.490	10.250	3.042	8.574	0.2187
555	44.270	33.800	9.829	3.281	8.606	0.2172
560	45.100	33.120	9.404	3.531	8.628	0.2156
565	45.920	32.450	8.979	3.791	8.638	0.2141
570	46.730	31.800	8.557	4.060	8.637	0.2125
575	47.520	31.160	8.137	4.339	8.626	0.2110
580	48.300	30.540	7.721	4.626	8.604	0.2095
585	49.050	29.940	7.309	4.920	8.573	0.2080
590	49.790	29.350	6.903	5.221	8.531	0.2065
595	50.500	28.790	6.502	5.526	8.480	0.2051
600	51.190	28.240	6.109	5.836	8.421	0.2037
605	51.860	27.710	5.724	6.149	8.354	0.2023
610	52.500	27.200	5.349	6.462	8.279	0.2009
615	53.120	26.720	4.983	6.776	8.198	0.1996
620	53.710	26.260	4.629	7.089	8.111	0.1983
625	54.270	25.820	4.287	7.398	8.020	0.1970
630	54.810	25.410	3.958	7.704	7.924	0.1959
635	55.310	25.020	3.643	8.003	7.825	0.1947
640	55.780	24.660	3.342	8.296	7.723	0.1936
645	56.220	24.330	3.057	8.582	7.620	0.1926
650	56.630	24.010	2.788	8.858	7.515	0.1916
655	57.010	23.730	2.534	9.124	7.410	0.1907
660	57.360	23.460	2.297	9.380	7.305	0.1898
665	57.680	23.230	2.075	9.625	7.201	0.1890
670	57.970	23.010	1.870	9.859	7.099	0.1883
675	58.230	22.820	1.681	10.080	6.997	0.1876
680	58.460	22.650	1.507	10.290	6.898	0.1870
685	58.670	22.510	1.349	10.490	6.800	0.1864
690	58.850	22.380	1.204	10.680	6.705	0.1859
695	59.010	22.270	1.073	10.850	6.613	0.1854
700	59.140	22.170	0.954	11.020	6.522	0.1850];
%T    H2(g)	H2O(g)	CH4(g)	CO(g)	CO2(g)	N2(g)
hscData=hscData(:,[1 2 4 5 6 3 7]);



%%%  Define the mixture model for reformate
% Use this vector order always H2, CH4, C0, CO2, H2O
fuel=MediumModel({'H2','CH4','CO','CO2','H2O','N2'});
% Assume Steam:Carbon=2.5
Z=[0 1 0 0 2.5 0.002]';
fuel.setZ(Z./sum(Z));
fuel.setT([hscData(:,1) ]+273.15);

nu=[    [3 -1 1  0 -1 0 ]' ...
        [1 0  -1 1 -1 0]'];
    
fuel.setNu(nu);
fuel.gibbs;
fuel.solveEq


figure
hold all
polyMeth = [-8.06168e-11,+3.04817e-7,-4.3751e-4,+0.296549,-81.903];
polyWgs = [1.70079e-11,-6.41497e-8,9.14923e-5,-6.07656e-2,16.576];
kMeth = 10.^(polyval(polyMeth,fuel.T));
kWgs = 10.^polyval(polyWgs,fuel.T);

plot(fuel.T-273.15,fuel.ln_kc(:,1))
plot(fuel.T-273.15,fuel.ln_kc(:,2))

plot(fuel.T-273.15,log(kMeth),':')
plot(fuel.T-273.15,log(kWgs),':')
legend('NASA K_c : Meth','NASA K_c : WGS','Leah K_c : Meth','Leah K_c : WGS')
title('Comparison of OLD Leah Polynomial fit and NASA Glenn')





% Comparison of ??HSc??Excel?? and NASA Glenn for Reforming

figure
hold all
ZDryEq=fuel.Zeq(:,[1:4 6])./repmat(sum(fuel.Zeq(:,[1:4 6]),2),1,5);
facCol=colormap(lines);
h1 = plot(fuel.T-273.15,100.*fuel.Zeq./repmat(sum(fuel.Zeq,2),1,6));
h2 = plot(hscData(:,1),hscData(:,2:7),'--');
set(h1(1),'Color',facCol(1,:));
set(h1(2),'Color',facCol(2,:));
set(h1(3),'Color',facCol(3,:));
set(h1(4),'Color',facCol(4,:));
set(h1(5),'Color',facCol(5,:));
set(h1(6),'Color',facCol(6,:));
set(h2(1),'Color',facCol(1,:));
set(h2(2),'Color',facCol(2,:));
set(h2(3),'Color',facCol(3,:));
set(h2(4),'Color',facCol(4,:));
set(h2(5),'Color',facCol(5,:));
set(h2(6),'Color',facCol(6,:));
legend({fuel.names{:} fuel.names{:}})
ylabel('Molar Composition [%]')
xlabel('Temperature [degC]')
title('Comparison of data in Excel(:) with NASA Glenn(-)')



% Check composition at 600degC
tReformateTest=600;
rH2Hsc=interp1(hscData(:,1),hscData(:,2),tReformateTest);
rH2MedMod=100*interp1(fuel.T,fuel.Zeq(:,1),tReformateTest+273.15);

assert(abs(rH2Hsc-rH2MedMod) < 1)

rCO2Hsc=interp1(hscData(:,1),hscData(:,5),tReformateTest);
rCO2MedMod=100*interp1(fuel.T,fuel.Zeq(:,4),tReformateTest+273.15);

assert(abs(rCO2Hsc-rCO2MedMod) < 1)

% Check composition at 550degC
tReformateTest=550;
rH2Hsc=interp1(hscData(:,1),hscData(:,2),tReformateTest);
rH2MedMod=100*interp1(fuel.T,fuel.Zeq(:,1),tReformateTest+273.15);

assert(abs(rH2Hsc-rH2MedMod) < 1)

rCO2Hsc=interp1(hscData(:,1),hscData(:,5),tReformateTest);
rCO2MedMod=100*interp1(fuel.T,fuel.Zeq(:,4),tReformateTest+273.15);

assert(abs(rCO2Hsc-rCO2MedMod) < 1)

%% Check equilibrium temperature calculation methods
tCatalyst=0*fuel.T+600+273.15;
[tEq,tATE]=equilibriumTemperature(fuel,fuel.Zeq,tCatalyst);

if max(abs(tEq-fuel.T*[1 1]))>1
    error('MediumModel:SelfTestReformate:ATECalc','Failure to correctly compute equilibrium temperature from a given composition')
end

if max(abs(tEq+tATE-tCatalyst*[1 1]))>1e-6
    error('MediumModel:SelfTestReformate:ATECalc','Failure to correctly compute ATE from a given composition')
end

%% check equilibrium method
 ln_k=fuel.equilibrium(fuel.Zeq,nu);

if max(abs((ln_k'-fuel.ln_kc)./fuel.ln_kc)) > 1e-3
    error('MediumModel:SelfTestReformate:lnkcalc','Computation of ln_k via 2 methods not self consistent')
end
disp('MediumModel.SelfTestReformate -- Test Passed')
