function tOut = findTFromH(me,h,swtPlot)
% Function to find T given h for a mixture of gases of compositions me.Z.
% Use 2 interp lookups to find the T given the value of h.
% h can be a column vector, in which case me.Z must be a single row or must 
% be the same row size as h.

% Store the original T
Z = me.Z;

% Size checks
h = h(:);
tOut = zeros(size(h));
if size(me.Z,1)>1 & length(h)==1
    h = repmat(h,1,size(me.Z,1));
    tOut = tOut';
elseif size(me.Z,1)>1 & length(h)>1
    h = h';
    tOut = tOut';
end

% Loop for each composition, a different h column for each one
for ctComp=1:size(me.Z,1)
    % Find outermost temperature range of gases data, if Z(gas)>0
    tMin = NaN; tMax = NaN;
    for ctGas=1:length(me.names)
        if me.Z(ctComp,ctGas)<eps
            continue;
        end
        tRange = cell2mat(me.gas.(me.names{ctGas}).tRange);
        tMinWc = max(tMin,min(tRange)); % inner, worst case gas ranges
        tMaxWc = min(tMax,max(tRange));
        tMin = min([tMin,tRange]);
        tMax = max([tMax,tRange]);
    end
    
    % Perform 2 vector interps to home in on answer.
    % This is more accurate than the previos fminsearch method & 2* quicker.
    T_V = linspace(tMin,tMax,200)';
    TWc_V = linspace(tMinWc,tMaxWc,200)';
    me_h = getH(me,T_V,Z(ctComp,:));
    me_hWc = getH(me,TWc_V,Z(ctComp,:));
    if min(h(:,ctComp))<me_hWc(1) | max(h(:,ctComp))>me_hWc(end)
        warning('MediumModel:findTFromH_OutOfRange','Beyond range of defined polynomial at least one gas - extrapolating.');
        tOut(:,ctComp) = interp1(me_hWc,TWc_V,h(:,ctComp),'linear','extrap');
    else
        if length(unique(me_h))<length(me_h)
            % Beyond the defined polynomial range for a gas, this can be expected
            tOut(:,ctComp) = NaN;
            continue;
        end
        tOut(:,ctComp) = interp1(me_h,T_V,h(:,ctComp),'linear','extrap');
        x1 = find(T_V<min(tOut(:,ctComp)),1,'last');
        x2 = find(T_V>max(tOut(:,ctComp)),1,'first');
        if isempty(x1)
            x1 = x2;
            x2 = x2+1;
        elseif isempty(x2)
            x2 = x1;
            x1 = x1-1;
        end
        T_V = linspace(T_V(x1),T_V(x2),200)';
        me_h = getH(me,T_V,Z(ctComp,:));
        tOut(:,ctComp) = interp1(me_h,T_V,h(:,ctComp),'linear','extrap');
    end
end
tOut = tOut(:);

% Optional plot - only for scalar case
if nargin>2 & swtPlot & size(h,1)==1
    T_V = (0:10:1500)'+273.15;
    h_V = getH(me,T_V,Z);
    figure;
    plot(h_V,T_V);grid on;hold on;
    plot(h(1),tOut,'r+','MarkerSize',12);
    xlabel('h [J/(kg.K)]');
    ylabel('T [K]');
end
end

function h = getH(me,T,Z)
% Function to evaluate only h - taken from props function
if size(Z,1)==1
    Z=repmat(Z,length(T),1);
end

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
