function SelfTestSimpleChp()
close all

oChp = CSimpleChp.Chp;
oChp.solve

oChp.testChp


%%
disp('MediumModel.SelfTestSimpleChp -- Test Passed')



end