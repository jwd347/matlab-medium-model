function solveEqQSearch(me)
% solveEqQSearch Solve for equilibrium condition
% See help MediumModel.setNu
%
% Uses reaction co-ordinate method
% Z = Zo+ facReacCoord*nu
% ln_Q =(ln Z)' * nu
% Aim: Find facReacCoord such that ln_Q=ln_k
% This can be achieved by finding the facReacCoord that minimises
% ln_k_err = ln_Q - ln_k
%          = ln(Zo+facReacCoord*nu)'*nu-ln_k
% Zeq is then given by
%   Zeq=Zo+ facReacCoord*nuEq
if isempty(me.ln_kc)|true % always do this
	me.gibbs;
end

me.Zeq=zeros(length(me.T),length(me.Z));
nu=me.nu;
strFminOpt=optimset(@fminsearch);
strFminOpt.TolFun =me.tol*1e-3 ;
strFminOpt.TolX=10^-6;
for ctT = 1:length(me.T)
	ln_k=me.ln_kc(ctT,:);
	if any(isnan(ln_k)) % added by ACM to prevent invalid Eq Concs
		disp('******* warning  - ln_k value not defined *************')
		me.Zeq(ctT,:) = ones(size(me.names))*NaN;
    else
            
		   facReacCoord0=zeros(1,size(nu,1)); % set rate of zero as initial guess. Has one entry per reaction
           hnd_k_err=@(facReacCoord) sum(abs(ComputeLnQFromCoord(me,facReacCoord,me.Z,nu)-ln_k)); % function to compute error between reaction quotient and equilibrium constant
		   [facReacCoordEq k_err]=fminsearch(hnd_k_err,facReacCoord0,strFminOpt);
		   me.Zeq(ctT,:)=ComputeZFromCoord(facReacCoordEq,me.Z,nu);
           ln_Q=ComputeLnQFromCoord(me,facReacCoordEq,me.Z,nu);
		  
		if abs(sum(me.Zeq(ctT,:))-1) > 1e-3
           warning('MediumModel:SolveEq','Composition vector does not sum to 1 (=%f)',sum(me.notCondensed(:,1).*me.Zeq(ctT,:)))
		end
		
	end
    
    % check whether a good solution was found:
	if max(abs(ln_Q-ln_k)./abs(ln_k)) > me.tol 
		warning('MediumModel:SolveEq','Equilibrium not found at %4.1f degC, point %d (%f>%f)',me.T(ctT)-273.15,ctT,max(abs(ln_Q(:)-ln_k(:))),me.tol);
        me.Zeq(ctT,:)=nan;
	end

end
me.moleToMassFractions
me.props
end


function lnQ=ComputeLnQFromCoord(me,facReacCoord,Z0,nu)

Z=ComputeZFromCoord(facReacCoord,Z0,nu);

% compute activity. Note that Z/sum(Z.*me.notCondensed) gives mole fractions of
% gaseous components
a=me.P*Z/(sum(Z.*me.notCondensed)* me.P0);

lnQ=log(max(a,0*a+1e-9))*(nu.*repmat(me.notCondensed,size(nu,1),1))'; % function to compute log reaction quotient with div by zero protection
 
end

function Z=ComputeZFromCoord(facReacCoord,Z0,nu)
% compute mole fraction resulting from facReacCoord and initial mole
% fraction Z0
Z= max(0*Z0,Z0+facReacCoord*nu)./sum(max(0,Z0+facReacCoord*nu)); 

if min(Z0+facReacCoord*nu)<0
    Z=nan(size(Z));
end

end