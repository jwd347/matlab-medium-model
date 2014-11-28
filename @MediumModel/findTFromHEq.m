function tOut = findTFromHEq(me,h,swtPlot)
% tOut = findTFromHEq(me,h,swtPlot)
% Function to find T given h for a mixture of gases of compositions me.Z after equilibrating 
% according to the reaction defined in me.nu.
% Use fminsearch to find the T given the value of h using solveEq at each step.
% NB this function takes much longer than findTFromH.
% me.T, me.Z and h should be the same row size.

if isempty(me.nu), error('MediumModel:findTFromHEq_NoNu','Error - no Medium nu set.'); end;

% Create strFminOpt to make the computation in solveEq faster.
strFminOpt = optimset(@fminsearch);
strFminOpt.TolFun=inf ;

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
    hnd= @(T) abs(h(ctT) - getHEq(me,T,Z(ctT,:),strFminOpt));
    [tOut(ctT),FVAL,EXITFLAG,OUTPUT] = fminsearch(hnd,tOutGuess);
end
    
% Optional plot - only for scalar case
if nargin>2 & swtPlot & length(h)==1
    T_V = (0:10:1500)'+273.15;
    h_V = getHEq(me,T_V,Z,strFminOpt);
    figure;
    plot(h_V,T_V);grid on;hold on;
    plot(h,tOut,'r+','MarkerSize',12);
    xlabel('h [J/(kg.K)]');
    ylabel('T [K]');
end

% Put back original T & Z
me.setTandZ(T,Z);

end

function h = getHEq(me,T,Z,strFminOpt)
% Function to evaluate h after equilibrating

me.setTandZ(T,Z);
me.solveEq([],strFminOpt);
h = me.h;

end
