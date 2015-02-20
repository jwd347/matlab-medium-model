function tOut = findTFromHEq(me,h,TAte,TLim,swtPlot)
% tOut = findTFromHEq(me,h,TAte,TLim,swtPlot)
% Function to find T given h for a mixture of gases of initial compositions me.Z after  
% equilibrating at the medium temperature-TAte according to the reaction defined in me.nu. 
% TLim (optional) limits the search between 2 values of T - defaults to [-inf inf]. 
% Uses fminsearch to find the T given the value of specific enthalpy h using solveEq at each step.
% NB this function takes much longer than findTFromH.
% me.T, me.Z and h should be the same row size. swtPlot only works for the scalar case.

if isempty(me.nu), error('MediumModel:findTFromHEq_NoNu','Error - no Medium nu set.'); end;

% Create strFminOpt to make the computation in solveEq faster.
strFminOpt = optimset(@fminsearch);
strFminOpt.TolFun=inf ;

% Set default TAte
if nargin<3 || isempty(TAte)
    TAte = 0;
end

% Set default T range
if nargin<4 || isempty(TLim)
    TLim = [-inf inf];
end

T=me.T;
if size(me.Z,1)==1
    Z=repmat(me.Z,size(me.T,1),1);
else
    Z=me.Z;
end

ctNoNan = find(~isnan(me.T));
tOut = T;
tOut(~ctNoNan) = 0;

for ctT=ctNoNan'
    tOutGuess = T(ctT);
    hnd= @(T) abs(h(ctT) - getHEq(me,T,Z(ctT,:),strFminOpt,TLim,TAte));
    [tOut(ctT),FVAL,EXITFLAG,OUTPUT] = fminsearch(hnd,tOutGuess);
end
    
% Optional plot - only for scalar case
if nargin>4 && swtPlot && length(h)==1
    T_V = (0:10:2000)'+273.15;
    h_V = getHEq(me,T_V,Z,strFminOpt,[-inf inf],TAte);
    figure;
    plot(h_V,T_V);grid on;hold on;
    plot(h,tOut,'r+','MarkerSize',12);
    xlabel('h [J/mol]');
    ylabel('T [K]');
end

% Put back original T & Z
me.setTandZ(T,Z);

end

function h = getHEq(me,T,Z,strFminOpt,TLim,TAte)
% Function to evaluate h after equilibrating

if T(1)<TLim(1) || T(1)>TLim(2)
    h = inf;
else
    me.setTandZ(T,Z);
    me.solveEq([],strFminOpt,TAte);
    h = me.h;
end

end
