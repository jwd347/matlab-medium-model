%% Introduction to the use of MediumModel class on Matlab
% The MediumModel class is used for finding the equilibrium state of a 
% mixture of species over a range of temperatures. The model uses data
% compiled by NASA in the paper:                                    
%   NASA Glenn Coefficients for Calculating
%   Thermodynamic Properties of Individual Species
%   
%   Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon
%   Glenn Research Center, Cleveland, Ohio
%  http://www.grc.nasa.gov/WWW/CEAWeb/TP-2002-211556.pdf
%

%% Instantiating a MediumModel Object

Example1 = MediumModel({'H2O'})
%%
% In this example, the model has only a single species. The model uses the 
% strMaster database, stored in the IdealGases file to find the relevant
% values using polynomial approximations over the required temperature range. 
%
% *Important values in the class* 
%
% Defined by the user:
%%
% * T - Temperature Range
% * Z - Initial molar fractions of each species   _Default-equal molar proportions_
% * X - Initial mass fractions of each species. _This is calculated if Z is
% specified and vice versa_
% * P0 -Atmospheric Pressure _Default - 1Bar_
% * P  -Pressure of reaction _Default - 1Bar_
% * nu - Stochiometry of chemical equations
% 
% 
% Calculated based on user inputs: (All of these are matrices of values
% over the temperature range)
%%
% * cp / cp_V - Specific heat capacity
% * h / h_V - Specific enthalpy
% * s / s_V - Specific entropy
% * mu / mu_V - Chemical potential
% * mm / mm_V - Molar mass
% * Zeq - Molar fractions of each species during the reaction
%
%  The _V element of the parameters such as cp_V indicates that it is a vector.
%  It differs from cp by that cp is the specific heat capacity of the
%  mixture as a whole, while cp_V is a matrix showing the specific heat
%  capacity of each individual species.
%  Zeq is calculated after calling the SolveEq method, the others are properties
%  obtained from the IdealGases database polynomial approximations.
% 


%%  Modelling Reaction Systems
% This example is a model of the Haber process reaction used for making
% ammonia. The model will be used to find the proportions of each gas over
% the temperature range -50-500°C
%

Example2 = MediumModel({'N2', 'H2', 'NH3'});
%%
% This creates an object with the specified species included in the model,
% and must include all products and reactants

%% Specify reaction conditions
% 
% The class requires 4 parameters to calculate how any reaction is
% expected to progress: Initial composition, stochiometry, temperature and pressure.

%% Set composition
% The initial proportions of the species can be set by masses or moles.
% This example uses moles. (Z)

Example2.setZ ([0.25, 0.75, 0]);

%%
% The proportions are done in the same order as they are given in the class
% definition and must sum to 1.00. (This example has 25% Nitrogen, 75% Hydrogen, 0% Ammonia)

%% Set Stochiometry
% The stochiometry of the chemical reactions must be specified in the model.
% In the example of the Haber process, the relevant reaction is 
%%
%
% $$N_2 + 3H_2 \leftrightarrow 2NH_3$$
Example2.setNu ([-1 -3 2]');
%%
% Reactants are given negative numbers, as they are used up, and products
% are given a positive number. This is also done using the same order as they were
% defined in the class. Note transpose operator, as each reaction is defined
% in its own column rather than row.
% Models with multiple reactions taking place can also be modeled using a
% nu matrix with additional columns, defining additional reactions. (See
% example 3)

%% Set Temperature
% The model is often used to find the composition change over a range of
% temperatures, so the temperature is set as a range of values. 

Example2.setT((-50:10:500)+273.15) ;

%%
% Here, the range is set from *-50-500°C* , in steps of 10°C. It must also be
% adjusted for the Kelvin temperature scale by adding 273.15 throughout.

%% (Optional) Set Pressure
% The reaction pressure is set as a default of 100,000Pa 1Bar,
% but can be changed by changing the property P of the object. Here it is
% changed to 5 bar. 
Example2.P=500000 ;
%%
% The pressure of the environment can also be changed
% from its default of 1Bar. Here it is modified to a slightly more accurate
% value of atmospheric pressure
Example2.P0=101300;

%% Class Processes
% The reaction calculations are done by two methods.
Example2.gibbs	;
Example2.solveEq;
%%
% *Example2.gibbs* uses the Gibbs equation, $$\Delta g=h-Ts$ to calculate the lowest free
% energy state, and this can be plotted using the command: *Example2.gibbsPlot* ;
% This command is executed in the section below, also showing the graphs it
% produces
% *Example2.solveEq* solves the relevant equations to find the equilibrium proportions so the end results can be
% plotted


%% Output data and plotting
%
% The Matrix *Example2.Zeq* holds the molar compositions of each species in
% each column of the matrix respectively, increasing in temperature, ready
% to be plotted against temperature using the code below.

figure					

plot(Example2.T-273.15,Example2.Zeq);
legend(Example2.names);
grid on;
xlabel('T (°C)')
ylabel('molar fraction (mol/mol C)')

Example2.gibbsPlot


%% Reforming example
% This example will use the reforming reaction, which has two reactions
% taking place side by side whilst including the inert nitrogen.
%%
%
% $$CH_4 + H_2O \leftrightarrow 3H_2 + CO$$
%
% $$CO + H_2O \leftrightarrow H_2 + CO_2$$


Example3 = MediumModel({'H2','CH4','CO','CO2','H2O', 'N2'});
Z=[0 1 0 0 2.8 0.2]';
Example3.setZ(Z./sum(Z));
%%
% The object is defined with all required species, and the initial molar
% compositions are set by initially defining the proportions of each
% species, then dividing through by the sum to normalise the vector to sum to one as
% required.

Example3.setT([300:10:750]+273.15);
nu=   [ [3 -1 1  0 -1 0]' ...     
       [1 0  -1 1 -1 0]']    ;           
Example3.setNu(nu);
%%
% The nu matrix now has 2 columns, each with 6 rows, to represent the two
% reactions taking place. Remembering products +ve, reactants -ve

Example3.gibbs;
Example3.solveEq;

plot(Example3.T-273.15,Example3.Zeq);
legend(Example3.names);
grid on;
xlabel('T (°C)')
ylabel('molar fraction (mol/mol C)')














