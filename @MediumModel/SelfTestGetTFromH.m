function SelfTestGetTFromH()
% Execute me.findTFromH  & me.findTFromHEq with combinations of vectors for
% me.T, me.Z and h.


%% findTFromH Tests
%  ================
oMed = MediumModel({'O2','N2'});

% oMed.setTandZ([400,500]',[0.5 0.5;0.2 0.8]);
oMed.setTandZ(400,[0.5 0.5]);
oMed.T
oMed.Z
T = oMed.findTFromH(3500,1)
%%
tic;
for k=1:100
    T = oMed.findTFromH([3500,4000]');
end
toc
Tans = T

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


%% findTFromHEq Tests
%  ==================
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
T = oMed.findTFromHEq(-100000,1)

%%
tic;
for k=1:2
    T = oMed.findTFromHEq(-100000);
end
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


