function props(me)
%             me.CheckShape;

swtDispPoly=0;

P=me.P;
T=me.T;
R=me.R;

%   C T R = a1T–2 + a2T–1 + a3 + a4T + a5T 2 + a6T3+ a7T 4 (1)
X_cp_M =  [T.^-2 T.^-1 ones(size(T)) T T.^2 T.^3 T.^4];
% Ho(T)/RT = –a1T –2 + a2lnT/T + a3 + a4T/2 + a5T 2/3 + a6T 3/4 + a7T 4/5 + b1/T (2)
X_h0_M =  [-T.^-2 log(T)./T ones(size(T)) T./2 (T.^2)./3 (T.^3)./4 (T.^4)./5 T.^-1];
%             ((-alow[1] + T*(blow[1] + alow[2]*Math.log(T) + T*(1.*alow[3] + T*(0.5*alow[4] + T*(1/3*alow[5] + T*(0.25*alow[6] + 0.2*alow[7]*T))))))/T)
% So(T)/R = –a1T–2/2 – a2T–1 + a3lnT + a4T + a5T 2/2 + a6T 3/3
% + a7T 4/4 + b2 (3)#
X_s0_M =  [-(T.^-2)./2 -T.^-1 log(T) T  (T.^2)./2 (T.^3)./3 (T.^4)./4 ones(size(T))];
cp_V=ones(length(T),length(me.names));
h_V=ones(length(T),length(me.names));
s_V=ones(length(T),length(me.names));
me.mm_V=zeros(1,length(me.names));
for ctGas =1:length(me.names)
    
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
        
        
        if swtDispPoly & ~isempty(ind)
            DispPoly(me,GasName,ctRange)
        end
        
    end
end
me.cp_V=cp_V.*R; 
me.h_V=h_V.*R; 
me.s_V=s_V.*R; 


% if isempty(me.Zeq)
%     % #TODO - tidy up
%     disp('This should never be reached.');
%     Z=repmat(me.Z,length(me.T),1);
% else
    Z=me.Zeq;
% end
notCondensed=repmat(me.notCondensed,size(Z,1),1);
Znorm=Z./repmat(sum(Z.*notCondensed,2),1,size(Z,2));
swtIsNan=isnan(Znorm);
a_V=max((me.P/ me.P0)*Znorm.*notCondensed,~notCondensed); %compute activity
a_V(swtIsNan)=NaN;
me.mu_V=me.h_V-(me.s_V.*repmat(me.T,1,length(me.names)))+R*repmat(me.T,1,length(me.names)).*log(a_V); 

if any(abs(sum(Z,2)-1)>0.0001)
    warning('MediumModel:props:Composition vector does not sum to 1(min = %f,max= %f',min(sum(Z,2)),max(sum(Z,2)))
end
me.cp = sum(me.cp_V.*Z,2);
me.h = sum(me.h_V.*Z,2);
me.s = sum(me.s_V.*Z,2);
me.mu = sum(me.mu_V.*Z,2);
me.aeq = a_V;
% So(T)/R = –a1T–2/2 – a2T–1 + a3lnT + a4T + a5T 2/2 + a6T 3/3 + a7T 4/4 + b2 (3)
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
