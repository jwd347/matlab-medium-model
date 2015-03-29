function SelfTestFindTFromH()
% Execute me.findTFromH  & me.findTFromHEq with combinations of vectors for
% me.T, me.Z and h.


disp('Running SelfTestFindTFromH - findTFromH.');

%% findTFromH Tests
%  ================
oMed = MediumModel({'O2','N2'});

% oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
T1 = 400;
oMed.setTandZ(T1,[0.5 0.5]);
h1 = oMed.h;
T2 = oMed.findTFromH(h1,1);
assert(abs(T1-T2)<0.001,['findTFromH failed at ' num2str(oMed.T) 'K, ' mat2str(oMed.Z)]);

%%
tic;
ctCalc = 100;
h1 = [3500:100:4000]';
for k=1:ctCalc
    T = oMed.findTFromH(h1);
end
tDuration = toc;
Tans = T
disp(['Duration for ' num2str(ctCalc) ' findTFromH vector calcs: ' num2str(tDuration)]);

%%
% Vector tests
oMed.setTandZ(400,[0.5 0.5]);
T = oMed.findTFromH(3500)
oMed.setTandZ([400,500]',[0.5 0.5]);
T = oMed.findTFromH(3500)
oMed.setTandZ(400,[0.5 0.5;0.2 0.8]);
T = oMed.findTFromH(3500)
oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
T = oMed.findTFromH(3500)

oMed.setTandZ(400,[0.5 0.5]);
T = oMed.findTFromH([3500,4000]')
oMed.setTandZ([400,500]',[0.5 0.5]);
T = oMed.findTFromH([3500,4000]')
oMed.setTandZ(400,[0.5 0.5;0.2 0.8]);
T = oMed.findTFromH([3500,4000]')
oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
T = oMed.findTFromH([3500,4000]')

disp('Vector tests complete.');

%% 
% Test with liquid water
oMed = MediumModel({'CH4','H2O','H2ObLb'});
oMed.setTandZ(283.15,[0.5 0.25 0.25]);
oMed.findTFromH(oMed.h,1)
oMed.setT(599);
oMed.findTFromH(oMed.h,1)

Tin = [0 1 10 50 99 100 101 110]'+273.15;
oMed.setT(Tin);
h = oMed.h;
TReturn = oMed.findTFromH(h);
TError = TReturn - Tin
assert(max(abs(TReturn-Tin))<0.001,'findTFromH failed to return all matching temperatures.');

% A warning is expected from cases below 273.15K.
% A difference is expected here for cases below 273.15K.
% To use this function with gases including H2ObLb below 273.15, separate
% the problem into cases above & below 273.15.
disp('Warnings expected...');
Tin = [-10 -1]'+273.15;
oMed.setT(Tin);
h = oMed.h;
TReturn = oMed.findTFromH(h);

%% findTFromHEq Tests
%  ==================
disp('Running SelfTestFindTFromH - findTFromHEq.');


lstGasFuel={'H2','CH4','CO','CO2','H2O'};
nu=[[3 -1  1  0 -1]; ...
    [1  0 -1  1 -1]];
Z=[0 1 0 0 2.8]; Z=Z/sum(Z);
oMed = MediumModel(lstGasFuel);
% oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
oMed.setTandZ(400,Z);
oMed.setNu(nu);
oMed.T
oMed.Z
[T,ctIter] = oMed.findTFromHEq(-120000,50,[],1)

%%
tic;
lstGasFuel={'N2','O2','H2O','H2ObLb'};
nu=[0 0 -1 1];
Z=[0.75 0.15 0.1 0]; Z=Z/sum(Z);
oMed = MediumModel(lstGasFuel);
oMed.setTandZ(400,Z);
oMed.setNu(nu);
oMed.T
oMed.Z
[T,ctIter] = oMed.findTFromHEq(-22200,50,[],1)
toc

%%
lstGasFuel={'H2','CH4','CO','CO2','H2O'};
nu=[[3 -1  1  0 -1]; ...
    [1  0 -1  1 -1]];
Z=[0 1 0 0 2.8]; Z=Z/sum(Z);
oMed = MediumModel(lstGasFuel);
% oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
oMed.setTandZ(400,Z);
oMed.setNu(nu);
tic;
T = oMed.findTFromHEq(-100000,0,0.1);
toc
Tans = T
%%
% Vector tests
% me.T, me.Z and h should be the same row size.
Z=[0 1 0 0 2.8]; Z=Z/sum(Z);
Z = [Z;0 0.3 0 0 0.7];
h = [-1;-1.2]*1e5;
oMed.setTandZ([400,500]',Z);
T = oMed.findTFromHEq(h)

%%
% Test with liquid water
oMed = MediumModel({'CH4','H2O','H2ObLb'});

Tin = [0 1 10 50 99 100 101 110]'+273.15;
oMed.setTandZ(Tin,[0.5 0.25 0.25]);
oMed.setNu([0 -1 1]);
oMed.solveEq;
oMed.setTandZ(Tin,oMed.Zeq);
h = oMed.h;
tic;
[TReturn,ctIter] = oMed.findTFromHEq(h,0,0.001);
toc
ctIter
TError = TReturn - Tin
assert(max(abs(TReturn-Tin))<0.0001,'findTFromH failed to return all matching temperatures.');

%% Test beyond range of gas data
oMed.setTandZ(-10+273.15,[0.5 0 0.5]);
oMed.setNu([0 -1 1]);
h1 = -1.85e5;
disp('Warning expected...');
TReturn1 = oMed.findTFromHEq(h1,0,1,1)
h2 = -0.82e5;
TReturn2 = oMed.findTFromHEq(h2,0,1,1)

