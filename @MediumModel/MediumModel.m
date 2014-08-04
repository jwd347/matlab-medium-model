classdef MediumModel < handle
    %MEDIUMMODEL NASA Glen Polynomial Gas Properties.
    %
    %   NASA Glenn Coefficients for Calculating
    %   Thermodynamic Properties of Individual Species
    %
    %   Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon
    %   Glenn Research Center, Cleveland, Ohio
    %
    %   This implementation by Mark Selby, Ceres Power Ltd, 2011.
    %
    %   All published substances(~2000) are implemented, but for typical needs
    %   the following substance names are useful
    %       N2,O2,CH4,CO,CO2,H2,H2O
    %
    %   All properties are reported in SI units with Mol the internal unit
    %   of "amount of substance". Molar Mass of the mixture is maintained
    %   by all methods so conversions to mass are conviently available.
    %
    %   Key methods (note most properties are protected and can only be set
    %   through methods, this is because shape is important to alot of the
    %   calcs).
    %   medium.setT --> sets the temeperature (vectors or scaler) in K
    %   medium.setZ --> sets the molar concentration(resets Zeq)
    %   medium.setX --> sets the mass comcentration
    %   medium.setNu--> sets the Stoichiometry Matrix for reacting fluids
    %   props   -->  calculate fluid properties based on Zeq
    %   plot    -->  plot the  fluid properties based on Zeq
    %   gibbs   -->  calculates the gibbs free energy for the specified
    %                conditions
    %   plotGibbs->  plots the above
    %   solveEq  ->  using the equilibrium constant calc'd from G_o, solve
    %                the mixture equilibrium at each temperature.
    %
    %   PRoperties
    %   h   --> specific enthalpy J/Mol
    %   s   --> specific entropy J/mol.K
    %   cp  --> specific heat capacity J/mol.K
    %   mm  --> molar mass of the mix g/mol
    %
    %   Example
    %
    %       medium=MediumModel({'N2','O2'}) % create an instance of a medium model
    %       medium.setZ([0.8 0.2]); % set the composition
    %       %    set the temperatures for which the properties are desired
    %       medium.setT(273.15:873.15);
    %       medium.plot % create a figure with the main properties displayed.
    %       % Calculate the change in Enthalpy between 2 temperatures
    %       tLow=273;
    %       tHigh=973;
    %       medium.setT([tLow tHigh]);
    %       QHeatMedium= diff(medium.h)/medium.mm % change in Enthalpy [J/kg/K]
    %       %   calculate the Specific Heat in the units of J/kg.K at the
    %       %   temperatures specified
    %       cp = medium.cp/medium.mm
    %
    %   % For more detailed examples including reformation see...
    %   [SVN]\OperatingGuides\MatlabTraining\MediumModel_Example.m
    %
    %   Equlibrium Example
    %       haber=MediumModel({'N2','H2','NH3'});
    %       haber.setT(([300 400:50:600])+273.15);
    %       haber.setZ([1 1 1]./3);
    %       haber.setNu([-1 -3 2]'); % see help for this function
    %       haber.gibbs;
    %       kc = exp(haber.ln_kc) ;
    %
    %       % Published in
    %       % Chemistry the Central Science" Ninth Ed., by: Brown, Lemay, Bursten, 2003, ISBN 0-13-038168-3
    %       kcVal = [4.34e-3 1.64e-4 4.51e-5 1.45e-5 5.38e-6 2.25e-6];
    %       figure
    %       semilogy(haber.T-273.15,kcVal,'o-',haber.T-273.15,kc,':.')
    %       legend('Published Value','Predicted Value')
    %       ylabel('K_c')
    %       xlabel('Temperature [degC]')
    %       title('Comparison of Equilibrium constants for the Haber process')
    %
    %   Water Example
    %   % create a model with water in all three states
    %   water=MediumModel({'H2Obsb','H2ObLb','H2O'});
    %   % SEE ALSO XSteam for full IF97 steam tables.
    %
    %   %% Reforming Example
    %%%  Define the mixture model for reformate
    % % Use this vector order always H2, CH4, C0, CO2, H2O
    % fuel=MediumModel({'H2','CH4','CO','CO2','H2O','N2'});
    % % Assume Steam:Carbon=2.5
    % Z=[0 1 0 0 2.8 0.002]';
    % fuel.setZ(Z./sum(Z));
    % nu=[    [3 -1 1  0 -1 0 ]' ...
    %         [1 0  -1 1 -1 0]'];
    %  fuel.setT((500:10:700)+273.15)
    % fuel.setNu(nu);
    % fuel.gibbs
    % fuel.gibbsPlot
    %
    %
    % fuel.solveEq
    % figure
    % subplot(2,1,1)
    % hold all
    % ZDryEq=fuel.Zeq(:,[1:4 6])./repmat(sum(fuel.Zeq(:,[1:4 6]),2),1,5);
    % facCol=colormap(lines);
    % h1 = plot(fuel.T-273.15,100.*fuel.Zeq./repmat(sum(fuel.Zeq,2),1,6));
    % legend({fuel.names{:}})
    % ylabel('Molar Composition [%]')
    % xlabel('Temperature [degC]')
    % rConvCH4 = 100.*(1-fuel.Zeq(:,[2 ])./sum(fuel.Zeq(:,[2:4]),2));
    % subplot(2,1,2)
    % plot(fuel.T-273.15, rConvCH4)
    % ylabel('Methane Conversion [%]')
    % xlabel('Temperature [degC]')
    
    
    
    
    properties
        tol=10^-7;
        gas;
        names={};
        % Single Gas properties
        cp_V=[]; % J/molK
        h_V =[]; % J/mol
        s_V=[]; % J/molK
        mu_V=[];
        % Mix properties
        cp=[]; % J/molK
        h =[]; % J/mol
        s=[]; % J/molK
        mu=[]; % J/mol chemical potential (or molar gibbs free energy of formation)
        P=10^5; % bulk pressure
       
        
        
        R=8.314510;
        P0=10^5; % standard pressure
        mm_V =[];
        mm;
        index;
    end
    
    properties(SetAccess=protected)
        %% Reaction properties
        ln_kc=[];
        nu=[];
        Zeq=[];
        aeq=[];
        
        Z =[]; % Molar composition
        X =[]; % Mass composition
        T=(273.15+(0:800))'; % K
        notCondensed=[];
    end
    methods
        [Teq,ATE]=equilibriumTemperature(me,Z,Tcatalyst)
        setTandZ(me,T,Z)
        props(me)
        solveEq(me,facReacCoord0)
        solveEqODE(me)
        g_reaction=gibbs(me)
        gibbsPlot(me)
        CheckShape(me)
        setNu(me,nu)
        plot(me)
        moleToMassFractions(me)
        massToMoleFractions(me)
        setT(me,T)
        setZ(me,Z)
        setX(me,X)
        tOut=findTFromH(me,h)
    end
    methods
        function medium = MediumModel(cellGas,varargin)
            persistent strMaster
            if isempty(strMaster) %~exist('strMaster','var')
                load IdealGases
            end
            %             warning('MAS possible issue exists with calculation of S --> has an impact on G=H-T*S --> kc USE WITH CARE')
            if ~iscell(cellGas)
                cellGas={cellGas};
            end
            medium.names=cellGas;
            medium.notCondensed=nan(length(cellGas),1);
            for ctGas =1:length(cellGas)
                try
                    medium.gas.(cellGas{ctGas})=strMaster.(cellGas{ctGas});
                    medium.index.(cellGas{ctGas})=ctGas;
                    medium.notCondensed(ctGas)=~strMaster.(cellGas{ctGas}).swtCondensed;
                catch ME
                    error('MediumModel:Gas Not Suppoted','Gas name not recognised : %s',cellGas{ctGas})
                end
            end
            
            medium.Z=ones(ctGas,1)/ctGas;
            
            medium.CheckShape;
            medium.Zeq=repmat(medium.Z',length(medium.T),1);
            medium.props;
            medium.moleToMassFractions;
            %             medium.setZ(medium.Z);
        end
    end
    
    methods(Static)
        dZ_V = reactionRate(t,Z,ln_k,nu,notCondensed)
        ln_k=equilibrium(Z,nu,notCondensed)
        SelfTest()
        SelfTestReformate()
        SelfTestHaber()
        SelfTestBoudouard()
        SelfTestElectrodePotentials()
        SelfTestTwoPhaseWater()
        SelfTestMethaneCombustion()
        SelfTestSimpleChp()
        SelfTestProfile()
        Tutorial()
        UserGuide()
    end
end

