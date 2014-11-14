function tOut = findTFromH(me,h,swtPlot)
% Function to find T given h for a mixture of gases of compositions me.Z.
% Use fminsearch to find the T given the function for h.
% me.T, me.Z and h should be the same row size.

if size(me.Z,1)==1
    Z=repmat(me.Z,size(me.T,1),1);
else
    Z=me.Z;
end

ctNoNan = find(~isnan(me.T));
tOut = me.T;
tOut(~ctNoNan) = 0;

for ctT=ctNoNan'
    tOutGuess = me.T(ctT);
    hnd= @(T) abs(h(ctT) - getH(me,T,Z(ctT,:)));
    [tOut(ctT),FVAL,EXITFLAG,OUTPUT] = fminsearch(hnd,tOutGuess);
end

% Optional plot - only for scalar case
if nargin>2 & swtPlot & length(h)==1
    T_V = (0:10:1500)+273.15;
    h_V = getH(me,T_V,Z);
    figure;
    plot(h_V,T_V);grid on;hold on;
    plot(h,tOut,'r+','MarkerSize',12);
    xlabel('h [J/(kg.K)]');
    ylabel('T [K]');
end
end

function h = getH(me,T,Z)
% Function to evaluate only h - taken from props function

% Ho(T)/RT = –a1T –2 + a2lnT/T + a3 + a4T/2 + a5T 2/3 + a6T 3/4 + a7T 4/5 + b1/T (2)
X_h0_M =  [-T.^-2 log(T)./T ones(size(T)) T./2 (T.^2)./3 (T.^3)./4 (T.^4)./5 T.^-1];
h_V=ones(length(T),length(me.names));
for ctGas =1:length(me.names)
    GasName = me.names{ctGas};
    ctTRanges = length(me.gas.(GasName).tRange);
    for ctRange = 1:ctTRanges
        tRange=me.gas.(GasName).tRange{ctRange};
        ind= find((tRange(1)<=T)&(T<=tRange(2))); 
        if isempty(ind)
            continue;
        end
        h_V(ind,ctGas) = (X_h0_M(ind,:)*[me.gas.(GasName).a{ctRange} me.gas.(GasName).b{ctRange}(1)]').*T(ind);
    end
end
h_V=h_V.*me.R; 
h = sum(h_V.*Z,2);
end

