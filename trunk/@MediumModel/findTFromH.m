function tOut = findTFromH(me,h,swtPlot)
% Function to find T given h for a mixture of gases.
% Use fminsearch to find the T given the function for h.

tOutGuess = 500;
hnd= @(T) abs(h - getH(me,T));
[tOut,FVAL,EXITFLAG,OUTPUT] = fminsearch(hnd,tOutGuess);

if nargin>2 && swtPlot
    T_V = [0:10:1500]+273.15;
    for ctT=1:length(T_V)
        h_V(ctT) = getH(me,T_V(ctT));
    end
    figure;
    plot(T_V,h_V);grid on;hold on;
    plot(tOut,h,'r+','MarkerSize',12);
end
end

function h = getH(me,T)
% Function to evaluate only h - taken from props function

% Ho(T)/RT = –a1T –2 + a2lnT/T + a3 + a4T/2 + a5T 2/3 + a6T 3/4 + a7T 4/5 + b1/T (2)
X_h0_M =  [-T.^-2 log(T)./T ones(size(T)) T./2 (T.^2)./3 (T.^3)./4 (T.^4)./5 T.^-1];
h_V=ones(length(T),length(me.names));
me.mm_V=zeros(length(me.names),1);
for ctGas =1:length(me.names)
    GasName = me.names{ctGas};
    me.mm_V(ctGas) = me.gas.(GasName).mm;
    ctTRanges = length(me.gas.(GasName).tRange);
    for ctRange = 1:ctTRanges
        tRange=me.gas.(GasName).tRange{ctRange};
        ind= find((tRange(1)<=T)&(T<=tRange(2))); %#ok<PROP>
        if isempty(ind)
            continue;
        end
        h_V(ind,ctGas) = (X_h0_M(ind,:)*[me.gas.(GasName).a{ctRange} me.gas.(GasName).b{ctRange}(1)]').*T(ind);
    end
end
h_V=h_V.*me.R; %#ok<PROP>

if isempty(me.Zeq)
    Z=me.Z';
else
    Z=me.Zeq;
end

if any(abs(sum(Z,2)-1)>0.0001)
    warning('MediumModel:props:Composition vector does not sum to 1(min = %f,max= %f',min(sum(Z,2)),max(sum(Z,2)))
end
h = sum(h_V.*Z,2);
end

