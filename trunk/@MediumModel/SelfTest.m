function SelfTest()
 clc
MediumModel.SelfTestReformate()
MediumModel.SelfTestHaber()

try
    MediumModel.SelfTestBoudouard()
catch ME
    ME.getReport
    disp('MediumModel.SelfTestBoudouard -- Test Failed')
end

MediumModel.SelfTestElectrodePotentials()
MediumModel.SelfTestTwoPhaseWater()
MediumModel.SelfTestMethaneCombustion()

    