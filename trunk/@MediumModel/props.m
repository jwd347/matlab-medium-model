function props(me)
%             me.CheckShape;

swtDispPoly=0;

T=me.T;
R=me.R;
ctNames = length(me.names);
ctT = length(T);

%   C T R = a1T–2 + a2T–1 + a3 + a4T + a5T 2 + a6T3+ a7T 4 (1)
X_cp_M =  [T.^-2 T.^-1 ones(size(T)) T T.^2 T.^3 T.^4];
% Ho(T)/RT = –a1T –2 + a2lnT/T + a3 + a4T/2 + a5T 2/3 + a6T 3/4 + a7T 4/5 + b1/T (2)
X_h0_M =  [-T.^-2 log(T)./T ones(size(T)) T./2 (T.^2)./3 (T.^3)./4 (T.^4)./5 T.^-1];
%             ((-alow[1] + T*(blow[1] + alow[2]*Math.log(T) + T*(1.*alow[3] + T*(0.5*alow[4] + T*(1/3*alow[5] + T*(0.25*alow[6] + 0.2*alow[7]*T))))))/T)
% So(T)/R = –a1T–2/2 – a2T–1 + a3lnT + a4T + a5T 2/2 + a6T 3/3
% + a7T 4/4 + b2 (3)#
X_s0_M =  [-(T.^-2)./2 -T.^-1 log(T) T  (T.^2)./2 (T.^3)./3 (T.^4)./4 ones(size(T))];
X_k_M =  [ones(size(T)) sqrt(T) T];
cp_V=ones(ctT,ctNames);
h_V=ones(ctT,ctNames);
s_V=ones(ctT,ctNames);
me.mm_V=zeros(1,ctNames);
dynVisc_V = nan(ctT,ctNames);
k_V = nan(ctT,ctNames);
for ctGas =1:ctNames
    
    GasName = me.names{ctGas};
    me_gas_Name = me.gas.(GasName); % dynamic field names are time consuming - do only once
    me.mm_V(ctGas) = me_gas_Name.mm;
    ctTRanges = length(me_gas_Name.tRange);
    for ctRange = 1:ctTRanges
        tRange=me_gas_Name.tRange{ctRange};
        ind= find((tRange(1)<=T)&(T<=tRange(2))); 
        if isempty(ind)
            continue;
        end
        cp_V(ind,ctGas) = X_cp_M(ind,:)*me_gas_Name.a{ctRange}';
        h_V(ind,ctGas) = (X_h0_M(ind,:)*[me_gas_Name.a{ctRange} me_gas_Name.b{ctRange}(1)]').*T(ind);
        s_V(ind,ctGas) = X_s0_M(ind,:)*[me_gas_Name.a{ctRange} me_gas_Name.b{ctRange}(2)]';
        
        if swtDispPoly && ~isempty(ind)
            DispPoly(me,GasName,ctRange)
        end
        
    end
    if isfield(me_gas_Name,'c')
        if me.notCondensed(ctGas)
            dynVisc_V(:,ctGas) = me_gas_Name.c(1)*T.^(me_gas_Name.c(3)/2)./(T+me_gas_Name.c(2));
        else
            dynVisc_V(:,ctGas) = me_gas_Name.c(1)*exp(me_gas_Name.c(2)./(R*T));
        end
        k_V(:,ctGas) = X_k_M*me_gas_Name.d';
    end
    
end
me.cp_V=cp_V.*R; 
me.h_V=h_V.*R; 
me.s_V=s_V.*R;

me.dynVisc_V = dynVisc_V*1e-6; % Data is in uPa.s
me.k_V = k_V;

Z=me.Zeq;

notCondensed=ones(size(Z,1),1)*me.notCondensed;
Znorm=Z./(sum(Z.*notCondensed,2)*ones(1,size(Z,2)));
swtIsNan=isnan(Znorm);
a_V=max((me.P/ me.P0)*Znorm.*notCondensed,~notCondensed); %compute activity
a_V(swtIsNan)=NaN;
me.mu_V=me.h_V-(me.s_V.*(me.T*ones(1,ctNames)))+R*(me.T*ones(1,ctNames)).*log(a_V); 

if any(abs(sum(Z,2)-1)>0.0001)
    warning('MediumModel:props','Composition vector does not sum to 1(min = %f,max= %f',min(sum(Z,2)),max(sum(Z,2)));
end
me.cp = sum(me.cp_V.*Z,2);
me.h = sum(me.h_V.*Z,2);
me.s = sum(me.s_V.*Z,2);
me.mu = sum(me.mu_V.*Z,2);
me.aeq = a_V;
% So(T)/R = –a1T–2/2 – a2T–1 + a3lnT + a4T + a5T 2/2 + a6T 3/3 + a7T 4/4 + b2 (3)

% Dynamic viscosity & thermal conductivity of the gaseous mixture - for some gaseous species only
% Viscosity mix formula (Herning and Zipperer) from http://petrowiki.org/Gas_viscosity
% If mixture has a gaseous species with no data, result is NaN.
swtNotCondensed = logical(me.notCondensed);
if any(me.notCondensed) && all(me.dynVisc_V(1,swtNotCondensed))>0
    sqrt_mm = Z*sqrt(me.mm_V.*me.notCondensed)';
    sqrt_mm_M = ones(ctT,1)*sqrt(me.mm_V.*me.notCondensed);
    me.dynVisc = sum(me.dynVisc_V(:,swtNotCondensed).*Z(:,swtNotCondensed).*sqrt_mm_M(:,swtNotCondensed),2)./sqrt_mm;
    % The Wilke Mixture Rule is rather involved - simplify to molar proportions
    me.k = sum(me.k_V(:,swtNotCondensed).*Znorm(:,swtNotCondensed),2);
else
    % All species condensed or some gaseous ones have no dynVisc & ThermCond data
    me.dynVisc = NaN;
    me.k = NaN;
end

% Density (ideal gas) for gaseous species only
mm_gas = Znorm*(me.mm_V.*me.notCondensed)';
if any(me.notCondensed)
   me.rho = me.P*mm_gas*0.001./(me.R.*me.T);
else
    me.rho = NaN;
end

% Evaluate mm & X
me.moleToMassFractions;

% Prandl number
me.pr = me.cp.*me.dynVisc./(me.mm.*me.k)*1e3;

end


function DispPoly(me,GasName,ctRange)
txtCoeff=num2str(me.gas.(GasName).a{ctRange}');
lstCpTerms={'T^-2', 'T^-1', '1', 'T', 'T^2', 'T^3', 'T^4'};
txt=[GasName ': Cp(T)='];
for i=1:length(lstCpTerms)
    txt=[txt '+' txtCoeff(i,:) '*' lstCpTerms{i}];
end
disp(strrep(txt,' ',''))

txtCoeff=num2str([me.gas.(GasName).a{ctRange} me.gas.(GasName).b{ctRange}(1)]');
lstHTerms={'-T^-2', 'log(T)/T' ,'1', 'T/2', '(T^2)/3', '(T^3)/4', '(T^4)/5' ,'T^-1'};
txt=[GasName ': H(T)='];
for i=1:length(lstHTerms)
    txt=[txt '+' txtCoeff(i,:) '*' lstHTerms{i}];
end
disp(strrep(txt,' ',''))
txtCoeff=num2str([me.gas.(GasName).a{ctRange} me.gas.(GasName).b{ctRange}(2)]');
lstSTerms={ '-(T^-2)/2', '-T^-1', 'log(T)', 'T',  '(T^2)/2', '(T^3)/3', '(T^4)/4', '1'};
txt=[GasName ': S(T)='];
for i=1:length(lstSTerms)
    txt=[txt '+' txtCoeff(i,:) '*' lstSTerms{i}];
end
disp(strrep(txt,' ',''))
end
