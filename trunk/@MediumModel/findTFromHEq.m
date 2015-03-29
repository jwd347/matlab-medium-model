function [tOut,ctIter] = findTFromHEq(me,h,tAte,hTol,swtPlot)
% tOut = findTFromHEq(me,h,tAte,hTol,swtPlot)
% Function to find T given h for a mixture of gases of initial compositions me.Z after  
% equilibrating at the medium temperature-tAte according to the reaction defined in me.nu. 
% tAte is an approach to equilibrium temperature.
% hTol (optional, default 1J/mol) is a tolerance for the search. 
% Uses an interpolation search to find the T given the value of specific enthalpy h using solveEq at each step.
% NB this function takes much longer than findTFromDH.
% me.T, me.Z and h should be the same row size. swtPlot only works for the scalar case.

if isempty(me.nu), error('MediumModel:findTFromHEq_NoNu','Error - no Medium nu set.'); end;

% Create strFminOpt to make the computation in solveEq faster.
strFminOpt = optimset(@fminsearch);
strFminOpt.TolFun=inf ;

% Set default tolerance
if nargin<4 || isempty(hTol)
    hTol = 1;
end

% Set default tAte
if nargin<3 || isempty(tAte)
    tAte = 0;
end

% Save starting values
T=me.T;
Zorig=me.Z;
if size(me.Z,1)==1
    Z=repmat(me.Z,size(me.T,1),1);
else
    Z=me.Z;
end

ctNoNan = find(~isnan(me.T));
tOut = T;
tOut(~ctNoNan) = 0;

swtWarnBeyondRange = 0;
ctIterMax = 50;
ctIter = zeros(size(tOut));
% Loop for all rows
for ctT=ctNoNan'
    % Find temperature range of all gases data
    tMin = NaN; tMax = NaN;
    tRanges = [-inf;inf]*ones(1,length(me.names));
    for ctGas=1:length(me.names)
        if Z(ctT,ctGas)<eps && all(me.nu(:,ctGas)==0)
            continue
        end
        tRange = cell2mat(me.gas.(me.names{ctGas}).tRange)';
        tRanges(:,ctGas) = [min(tRange);max(tRange)];
        tMin = max(tMin,min(tRanges(:,ctGas))); % inner, narrowest gas ranges
        tMax = min(tMax,max(tRanges(:,ctGas)));
    end

    % Find h at min T
    me.setTandZ(tMin,Z(ctT,:));
    me.solveEq([],strFminOpt,tAte);
    h_V = me.h;
    t_V = tMin;
    if h(ctT)<=me.h && any(me.Zeq(tRanges(1,:)==tMin)>0.01)
        % Extrapolate beyond low end
        me.setTandZ(tMin+1,Z(ctT,:));
        me.solveEq([],strFminOpt,tAte);
        tOut(ctT) = interp1([h_V me.h],[tMin tMin+1],h(ctT),'linear','extrap');
        swtWarnBeyondRange = 1;
        continue
    end
    
%     hndF=figure;plot(h_V,t_V,'*');hold on; % debug plot

    % Find h at max T
    me.setTandZ(tMax,Z(ctT,:));
    me.solveEq([],strFminOpt,tAte);
    hNew = me.h;
    tNew = tMax;
    if h(ctT)>=me.h && any(me.Zeq(tRanges(2,:)==tMax)>0.01)
        % Extrapolate beyond high end
        me.setTandZ(tMax-1,Z(ctT,:));
        me.solveEq([],strFminOpt,tAte);
        tOut(ctT) = interp1([me.h hNew],[tMax-1 tMax],h(ctT),'linear','extrap');
        swtWarnBeyondRange = 1;
        continue
    end
    
    % If between limits, use pchip, else use linear which can extrapolate
    if h(ctT)<hNew && h(ctT)>h_V
        txtMethod = 'pchip';
    else
        txtMethod = 'linear';
    end
    
    % Search for tOut(ctT) at h(ctT)
    while abs(h(ctT)-hNew)>hTol
        
        [h_V,ind] = sort([h_V;hNew]);
%         figure(hndF);plot(hNew,tNew,'*'); % debug plot
        t_V = [t_V;tNew];
        t_V = t_V(ind);
        [h_U,indU] = unique(h_V);
        t_U = t_V(indU);
        tNew = interp1(h_U,t_U,h(ctT),txtMethod,'extrap');
        
        ctIter(ctT) = ctIter(ctT)+1;
        
        hNew = getHEq(me,tNew,Z(ctT,:),strFminOpt,tAte);
        
        if  ctIter(ctT)==ceil(ctIterMax/2)
            % Getting nowhere - maybe inaccurate points exist so clear out a bit
            % Replace all points within 0.1% of latest tNew with re-evaluated point at tNew
            warning('MediumModel:findTFromHEq_IterMid',['Iteration reached clear-out at ' int2str(ctIter(ctT)) ' iterations.']);
%             figure;plot(dH_V,t_V,'*'); hold on;% debug plot
%             title(['findTFromHEq_IterMid searching at h=' num2str(h(ctT))]);
%             vline(h(ctT),'r:');
%             plot(hNew,tNew,'ro');
            threshTV = (max(t_V)-min(t_V))/1000;
            indKeep = find(t_V<tNew-threshTV | t_V>tNew+threshTV);
            t_V = t_V(indKeep);
            h_V = h_V(indKeep);
        elseif  ctIter(ctT)>ctIterMax
            warning('MediumModel:findTFromHEq_IterMax','Iteration reached maximum.');
%             figure;plot(h_V,t_V,'*'); hold on;% debug plot
%             title(['findTFromHEq_IterMax searching at h=' num2str(h(ctT))]);
%             vline(h(ctT),'r:');
%             plot(hNew,tNew,'ro');
            break;
        end
        
    end

    tOut(ctT) = tNew;
end

if swtWarnBeyondRange
    warning('MediumModel:findTFromHEq_OutOfRange','Beyond range of defined polynomial at least one gas - extrapolating.');
end
    
% Optional plot - only for scalar case
if nargin>4 && swtPlot && length(h)==1
    % Add detail close to the result
    T_V = sort([linspace(tMin,tMax,50),...
        linspace(max(tOut-100,tMin),min(tOut+100,tMax),25),...
        linspace(max(tOut-10,tMin),min(tOut+10,tMax),25)]);
    h_V = getHEq(me,T_V,Z,strFminOpt,tAte);
    figure;
    plot(h_V,T_V,'-','MarkerSize',8);grid on;hold on;
    plot(h,tOut,'r+','MarkerSize',12);
    xlabel('h [J/mol]');
    ylabel('T [K]');
end

% Put back original T & Z
me.setTandZ(T,Zorig);

end

function h = getHEq(me,T,Z,strFminOpt,tAte)
% Function to evaluate h after equilibrating

me.setTandZ(T,Z);
me.solveEq([],strFminOpt,tAte);
h = me.h;

end
