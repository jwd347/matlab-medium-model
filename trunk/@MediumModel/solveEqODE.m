function solveEq(me)
% SOLVEEQ Solve for equilibrium condition
% See help MediumModel.setNu
if isempty(me.ln_kc)|true % always do this
    me.gibbs;
end

me.Zeq=zeros(length(me.T),length(me.Z));
nu=me.nu;
notSolids = me.notSolids;

for ctT = 1:length(me.T)
    ln_k=me.ln_kc(ctT,:);
    if any(isnan(ln_k)) % added by ACM to prevent invalid Eq Concs
        disp('******* warning  - ln_k value not defined *************')
        me.Zeq(ctT,:) = ones(size(me.names))*NaN;
    else
        
        hndReac = @(t,Z) me.reactionRate(t,Z,ln_k,nu,notSolids);
        optionsSolve = odeset('RelTol',1e-9,'AbsTol',1e-9);
        [tiVec,Y]=ode15s(hndReac,[0 4000000],me.Z,optionsSolve);
        me.Zeq(ctT,:)=Y(end,:)./sum(Y(end,:));
        if abs(sum(me.Zeq(ctT,:))-1) > 1e-3
            warning('MediumModel:SolveEq','Composition vector does not sum to 1 (=%f)',sum(me.Zeq(ctT,:)))
        end
        % #FIXME Needs to take account of "notSolids" to avoid throwing
        % this warning
        ln_k_val2=log(prod(repmat(me.Zeq(ctT,:),size(me.nu',1),1).^(me.nu'),2));
        ln_k_val=me.equilibrium(me.Zeq(ctT,:),me.nu,notSolids);
        if ~nearly(ln_k_val2,ln_k_val,0.00001)
            warning('MediumModel:SolveEq','Composition Solution does not have correct equilibrium position')
            disp('Theory, Actual')
            disp(num2str([ln_k_val2,ln_k_val],'\t%+8.6g'))
        end
    end
    
    if max(abs(ln_k_val-ln_k')) > me.tol
        try
            warning('MediumModel:SolveEq','Equilibrium not found at %4.1f degC, point %d (%f>%f)',me.T(ctT)-273.15,ctT,max(abs(ln_k_val-ln_k')),me.tol);
        end
        %                     figure
        %                     plot(tiVec,Y)
        %                     legend(me.names)
    end
    
end
me.moleToMassFractions
me.props
end
