function solveEq(me,facReacCoord0)
% SOLVEEQ Solve for equilibrium condition
% See help MediumModel.setNu
%
% Define reaction co-ordinate to be facReacCoord in the following
% N= N0+facReacCoord*nu; 
% where nu is the stoichiometry matrix, N is the moles of each species and No is the initial number of moles.
% The approach taken here is to search over all facReacCoord to find the
% final composition that minimises the Gibbs free energy, subject to the constraint
% that all elements of N must non non-negative. This is the
% equilibrium position because any change either side of this point would
% lead to an increase in Gibbs energy and therefore could not occur
% spontaneously.


% Process optional initial guess for the reaction co-ordinate
if nargin==1
    facReacCoord0=zeros(1,size(me.nu,1)); % by default set reaction co-ordinate to zero as initial guess. Has one entry per reaction
else
    facReacCoord0=reshape(facReacCoord0(:),1,size(me.nu,1));
end

% initialise equilibiurm concentration at zero
me.Zeq=zeros(length(me.T),size(me.Z,2));

% alias nu matrix
nu=me.nu;

% setup options for solver 
strFminOpt=optimset(@fminsearch);
strFminOpt.TolFun =inf ;
strFminOpt.TolX=me.tol;

% iterate over each temperature
for ctT = 1:length(me.T)
    if size(me.Z,1)>1
        % Range of starting compositions, one for each temperature
        Z = me.Z(ctT,:);
    else
        % Single starting composition for all temperatures
        Z = me.Z;
    end
    % Define anonymous function that computes gibbs energy as a function of reaction coord 
    hndGibbsFromReacCoord=@(facReacCoord)...
        ComputeGibbs(me,ComputeNFromReacCoord(Z,facReacCoord,nu),...
        me.h_V(ctT,:),me.s_V(ctT,:),me.T(ctT,:)) ;
    % find the reaction co-ord that minimises gibbs energy (i.e. the
    % equilibrium reaction co-ordinate)
    [facReacCoordEq Gmin ctFminExitFlag]=fminsearch(hndGibbsFromReacCoord,facReacCoord0,strFminOpt);
    
    % compute equilibrium composition from equilibrium reaction co-ord
    me.Zeq(ctT,:)=ComputeZFromCoord(facReacCoordEq,Z,nu)';
    
    % warn if solution not found
    if ctFminExitFlag<1
        warning('MediumModel:SolveEq','Equilibrium not found at %4.1f degC, point %d ',me.T(ctT)-273.15,ctT)
        me.Zeq(ctT,:)=nan;
    end
    
    % warn if composition vector doesn't add to 1
    if abs(sum(me.Zeq(ctT,:))-1) > 1e-3
        warning('MediumModel:SolveEq','Composition vector does not sum to 1 (=%f)',sum(me.notCondensed(:,1).*me.Zeq(ctT,:)))
    end
    
    
    
end
me.moleToMassFractions
% compute properties for equilibrium concentrations
me.props
end

function G=ComputeGibbs(me,N,h_V,s_V,T)
% compute gibbs energy as a function of number of moles of species (N) , 
% the enthalpy vector (h_V) and the entropy vector (s_V) and the
% temperature


% compute mole fractions of gaseous species from number of moles
ZnotCond=N./sum(N.*me.notCondensed);
swtIsNan=isnan(ZnotCond);

%compute activity
a_V=max((me.P/ me.P0)*ZnotCond.*me.notCondensed,~me.notCondensed); 
a_V(swtIsNan)=NaN;

% compute chemical potential
mu_V=h_V-(s_V.*T)+me.R*T*log(a_V); %#ok<PROP>
swtIsFinite=isfinite(mu_V);

% remove non-finite elements of mu_V resulting from zero activities and
% compute the total gibbs energy
G=mu_V(swtIsFinite)*N(swtIsFinite)';

if any(isnan(N))
  %  G=nan;
end
end

function N=ComputeNFromReacCoord(N0,facReacCoord,nu)
% compute number of moles of products fraction resulting from facReacCoord
% and initial moles N0 
N= N0+facReacCoord*nu; 

% return NaN if any element has become negative
if min(N)<0
    N=nan(size(N0));
end


end

function Z=ComputeZFromCoord(facReacCoord,Z0,nu)
% compute mole fraction resulting from facReacCoord and initial mole
% fraction Z0
Z= max(0*Z0,Z0+facReacCoord*nu)./sum(max(0,Z0+facReacCoord*nu)); 

if min(Z0+facReacCoord*nu)<0
    Z=nan(size(Z));
end

end