%% Introduction to the use of MediumModel class on Matlab
% The MediumModel class is used for finding the equilibrium state of a 
% mixture of species over a range of temperatures

%%
%% Class definition
Example2 =MediumModel({'N2', 'H2', 'NH3'});
%%
% This creates a class with the specified species included in the model,
% including all potential products and reactants


%% 
%% Specify reaction conditions
% 
%
% The class requires 4 parameters to calculate how any reaction is
% expected to progress: Initial composition, stochiometry, temperature and pressure.

%% Set composition
% The initial proportions of the species can be set by masses or moles

Example2.setZ ([0.25, 0.75, 0]);

%%
% The proportions are done in the same order as they are given in the class
% definition and must sum to 1.00. (25% Nitrogen, 75% Hydrogen)

%% Set Stochiometry
% The stochiometry of the chemical reactions must be specified in the model.
% In the example of the Haber process, the relevant reaction is 
% *N2 + 3 H2 < = > 2 NH3* .
Example2.setNu ([-1 -3 2]');
%%
% Reactants are given negative numbers, as they are used up, and products
% positive number. This is also done using the same order as they were
% defined in the class. Note transpose operator as each reaction is defined
% it its own column rather than row.
% Models with multiple reactions taking place can also be modeled using a
% nu matrix with additional columns, defining additional reactions

%% Set Temperature
% The model is often used to find the composition change over a range of
% temperatures, so the temperature is set as a range of values. 

Example2.setT((-50:10:300)+273.15) ;

%%
% Here, the range is set from *-100-300°C* , in steps of 10°C. It must also be
% adjusted for the Kelvin temperature scale by adding 273.15 throughout.

%% (Optional) Set Pressure
% The reaction pressure is set as a default of 100,000Pa (atmospheric pressure),
% but can be changed by changing the property P0 of the object.

Example2.P0=20000000 ;

%% Class Processes
% The reaction calculations are done by two methods.
Example2.gibbs	;
Example2.solveEq;
%%
% *Example2.gibbs* uses the Gibbs equation, G=H-TS to calculate lowest free
% energy state, and this can be plotted using the command: *Example2.gibbsPlot* ;
% This command is executed in the section below, also showing the graphs it
% produces
% *Example2.solveEq* solves the relevant equations so the end results can be
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
shg














