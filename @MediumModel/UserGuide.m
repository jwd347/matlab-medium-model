%%  MediumModel Class - User Guide
% The MediumModel class is a tool for modelling and analysing mixtures of
% reacting ideal gases over a range of temperatures, pressures and
% compositions. It can be used to calculate various basic thermodynamic
% properties of individual gases and gas mixtures.
%
% Reaction system equations can be defined and solved to find equilibrium
% compositions for a reaction (or set of reactions) at a range of
% temperatures.
%
% The model uses data compiled by NASA in the following paper:
%
% *NASA Glenn Coefficients for Calculating Thermodynamic Properties of
% Individual Species.* _Bonnie J. McBride, Michael J. Zehe, and Sanford
% Gordon_ Glenn Research Center, Cleveland, Ohio
% http://www.grc.nasa.gov/WWW/CEAWeb/TP-2002-211556.pdf
%
% All published substances(~2000) are implemented, but for typical needs
% the following substance names are useful *N2,O2,CH4,CO,CO2,H2,H2O*
%
% All properties are reported in SI units with mol the internal unit of
% "amount of substance". Molar Mass of the mixture is maintained by all
% methods so conversions to mass are conveniently available.
  

%% Class Properties
% *User-defined:  _(Using standard SI units)_*
%%
% <html>
% <table border=1><tr><td><b>Property</b></td><td><b>Units</b></td><td><b>Description</b></td><td><b>Default</b></td></tr>
% <tr><td><b>T</b></td><td>K</td><td>The temperature range defined for the medium model</td><td>273.15K -> 1073.15K in 1K steps</td></tr>
% <tr><td><b>Z</b></td><td>mol/mol</td><td>Initial molar fractions of each species in the medium model</td><td>equal fractions for all species</td></tr>
% <tr><td><b>X</b></td><td>Kg/Kg</td><td>Initial mass fractions of each species. This is calculated if Z
% is specified and vice versa</td><td>equal fractions for all species</td></tr>
% <tr><td><b>nu</b></td><td>n/a</td><td>Stoichiometry of defined chemical equations</td><td>empty</td></tr>
% <tr><td><b>P0</b></td><td>Pa</td><td>Atmospheric Pressure</td><td>100,000</td></tr>
% <tr><td><b>P</b></td><td>Pa</td><td>Reaction Pressure</td><td>100,000</td></tr></table>
% </html>
%
%%
% *Calculated:*
%%
% <html>
% <table border=1><tr><td><b>Property</b></td><td><b>Units</b></td><td><b>Description</b></td><td><b>Default</b></td></tr>
% <tr><td><b>cp / cp_V</b></td><td>J/mol.K</td><td>Specific heat capacity</td><td>n/a</td></tr>
% <tr><td><b>h / h_V</b></td><td>J/mol</td><td>Specific enthalpy</td><td>n/a</td></tr>
% <tr><td><b>s / s_V</b></td><td>J/mol.K</td><td>Specific entropy</td><td>n/a</td></tr>
% <tr><td><b>mu / mu_V</b></td><td>J/mol</td><td>Chemical potential</td><td>n/a</td></tr>
% <tr><td><b>mm / mm_V</b></td><td>g/mol</td><td>Molar mass</td><td>n/a</td></tr>
% <tr><td><b>Zeq</b></td><td>mol/mol</td><td>Molar fractions of each species at equilibrium over the
% temperature range</td><td>n/a</td></tr></table>
% </html>
%
%%
% *NB* The model uses the strMaster database, stored in the IdealGases file to
% find the relevant values using polynomial approximations over the
% required temperature range.

%%
% The '_V' element of the parameters such as cp_V indicates that it is a
% vector. It differs from cp by that cp is the specific heat capacity of
% the mixture as a whole, while cp_V is a matrix showing the specific heat
% capacity of each individual species. Zeq is calculated after calling the
% SolveEq method, the others are properties obtained from the IdealGases
% database polynomial approximations.
% 
%% Class Methods
% Methods are used in the class to set values of properties, perform
% calculations and plot graphs:
%
% *Property Set Methods*
%
% * *setT*     _-sets the object temperature range_
% * *setP*     _-sets the reaction and atmospheric pressures_
% * *setNu*    _-sets the object's chemical reaction stoichiometry_
% * *setZ*     _-sets the object's composition by moles_
% * *setX*     _-sets the object's composition by mass_
%
% *Calculation Methods*
%
% * *gibbs* _- uses the Gibbs equation, g=h-Ts to calculate the enthalpy
% change of reaction for the specified stoichiometry so the lowest free
% energy state can be calculated._
% * *solveEq* _- solves the relevant Gibbs equations to find the lowest
% free energy state and hence equilibrium proportions over the temperature
% range so the end results can be plotted._
% * *props* _- used to access the IdealGases database and uses the
% polynomial approximations to find the values of cp, h and s for the
% gases._
% * *moleToMassFractions* _- calculates relevant value by mass, if the user
% specifies values by moles_
% * *massToMoleFractions*  _- calculates relevant value by moles, if the
% user specifies values by mass_
% * *findTFromH* _- calculates temperatures when specific enthalpies are specified_
%
%
% |gibbs| and |solveEq| are executed in Example2 and Example3 below. The
% remaining calculation methods are internal methods that are not directly
% called by the user. |moleToMassFractions, massToMoleFractions| and
% |findTFromH| convert from one property to another as described, depending
% on which properties were initially defined by the user.
%
%
% *Plotting Methods*
%
% * plot
% * gibbsplot
%
% |plot| plots 4 graphs in one window, plotting each gas's cp, h, s and
% composition against temperature over the temperature range. |gibbsplot|
% produces two graphs, plotting free energy and equilibrium constant
% against temperature. |gibbsplot| is used in example2. |plot| is used in
% example3.
%
%

%% Example 1 - Instantiating a Simple MediumModel Object
% In this example, the model contains only a single species (steam).
Example1 = MediumModel({'H2O'});
Example1.setT([100:100:500]'+273.15);
Example1

%%
% The MediumModel object which is instantiated contains the set of
% thermodynamic properties for the defined species , as extracted from the
% NASA thermodynamic properties database. The media is defined for a range
% of temperatures from 373.115K to 773.15K (100 to 500°C).
% Below is a plot of the variation in heat capacity over the 5 temperature
% values.
plot(Example1.T,Example1.cp,'--o');xlabel('T (°C)');ylabel('cp(J/mol K');set(gca,'Xtick',Example1.T); 
%%  Example 2 - Modelling a Simple Reaction System
% This example demonstrates a model of the Haber process reaction, used for
% making ammonia.
%
% $$N_2 + 3H_2 \leftrightarrow 2NH_3$$
%%
% The model will be used to find the equilibrium compositions of each gas
% during the reaction process over the temperature range -50 to 500°C.

Example2 = MediumModel({'N2', 'H2', 'NH3'});
%%
% This creates an object with the specified species included in the model,
% and must include all products and reactants.
%
% *Specify reaction conditions*
% 
% The MediumModel class requires 4 parameters (Initial composition, 
% stoichiometry, temperature and pressure)
% to be set in order to calculate how any reaction will
% to progress as follow: 
%
% *setZ* (Composition)
%
% The initial proportions of the species are set by mass or molar
% fractions.

Example2.setZ ([0.25, 0.75, 0]);

%%
% The proportions are defined in the same order as they are given in the
% class definition, as stored by the 'names' property and must sum to 1.00.
% (This example has 25% Nitrogen, 75% Hydrogen, 0% Ammonia)
%
% *setNu* (stoichiometry)
%
% The stoichiometry of the chemical reactions must be specified in the
% model.

Example2.setNu ([-1, -3, 2]);
%%
% Reactants are given negative numbers, as they are used up, and products
% are given a positive number. This is also done using the same order as
% they were defined in the class.  Models with multiple reactions taking
% place can also be modelled using a nu matrix with additional columns,
% defining additional reactions. (See Example 3)
%
% *setT* (Temperature)
%
% The object can be used to find the composition change during the reaction
% over a range of temperatures by setting temperature property of the
% model.

Example2.setT((-50:10:500)'+273.15) ;

%%
% Here, the range is set from -50 to 500°C , in steps of 10°C.
%
% *setP* (Pressure)
%
% The reaction and atmospheric pressures are set by default to 100,000Pa,
% but can be modified by changing the property P of the object.
% Here, the reaction pressure is changed to 500,000Pa.
Example2.setP(500000);

%%
% *Solve Equilibrium Conditions*
% 
% The following commands calculate the equilibrium conditions for the
% defined system, populating the Zeq (equilibrium composition) property.
Example2.gibbs;
Example2.solveEq;
%%
% *Output data and plotting*
%
% The property *Example2.Zeq* holds the molar compositions of each species
% in each column of the matrix respectively, against temperature (as
% defined in  *Example2.T*). The output compositions can be plotted against
% reaction temperature as shown below:

figure					
plot(Example2.T-273.15,Example2.Zeq);
legend(Example2.names);
grid on;
xlabel('T (°C)')
ylabel('molar fraction (mol/mol C)')

%% Example 3 - Steam Methane Reforming
% This example will model the steam methane reforming reaction, which has two
% reactions taking place simultaneously.
%
% $$CH_4 + H_2O \leftrightarrow 3H_2 + CO$
%
% $$CO + H_2O \leftrightarrow H_2 + CO_2$

Example3 = MediumModel({'H2','CH4','CO','CO2','H2O'});
Z=[0 1 0 0 2.8];
Z=Z./sum(Z);
Example3.setZ(Z);
%%
% The object is defined with all required species, and the initial molar
% compositions are set by  defining the proportions of each
% species, then dividing through by the sum to normalise the vector to sum
% to one as required. In this case, a steam to carbon ratio of 2.8 is
% defined at the reformer input.

Example3.setT([300:10:750]'+273.15);
nu=   [ [3 -1 1  0 -1]; ...
      [1 0  -1 1 -1]]    ;
Example3.setNu(nu);
%%
% The nu matrix now has 6 columns, each with 2 rows, to represent the two
% reactions taking place. Note: products +ve, reactants -ve. 
% Once again  the equilibrium composition of the reaction against a
% range of temperatures can be found using the |solveEq| method.
%%
% The plotting method |plot|, is used to create plots of cp, s, h and
% equilibrium composition against temperature:
%%
Example3.gibbs;
Example3.solveEq;
Example3.plot