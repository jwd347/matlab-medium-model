function setTandZ(me,T,Z)
me.T=T;
me.Z=Z;
me.CheckShape;
if size(Z,1)==1
    me.Zeq=repmat(me.Z,length(me.T),1);
else
    me.Zeq = me.Z;
end

me.moleToMassFractions;
me.props
end