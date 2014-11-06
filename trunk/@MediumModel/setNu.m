function setNu(me,nu)
% SETNU Set the stoichiomtery matrix
% nu is an r by s matrix describing the equilibrium reaction
% system containg r reactions on s species
% -ve numbers == left hand side
% +ve numbers == right hand side
% 0 == not involved in reaction
% % Example : Reforming
% CH4+ H2O <--> 3H2 +CO     (1)
% CO + H2O <--> H2 + CO2    (2)
% fuel=MediumModel({'H2','CH4','CO','CO2','H2O','N2'});
% nu=[    [3 -1 1  0 -1 0 ]; ...
%         [1 0  -1 1 -1 0]];
% fuel.setNu(nu);
% fuel.gibbs(nu);
% fuel.gibbsPlot();
% fuel.solveEq();
% plot(fuel.T-273.15,fuel.Zeq(:,1).*100)
% xlabel('Temperature [degC]')
% ylabel('H_2 Fraction [%]')
me.nu=nu;
me.CheckShape;
ctSize= size(me.nu);
if ctSize(2)~=length(me.names)
    error('Stoichiometry Matrix must have the same number of columns as species in the model')
end
me.props;

end