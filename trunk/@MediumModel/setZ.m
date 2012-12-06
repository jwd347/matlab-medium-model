function setZ(me,Z)
if length(Z)~=length(me.names)
    error('MediumModel:setZ','Must have correct number of species')
end
if abs(sum(Z)-1) > 1e-3
    warning('MediumModel:setZ','Composition vector does not sum to 1 (=%f)',sum(Z))
end

me.Z=Z;
me.CheckShape;
me.Zeq=repmat(me.Z',length(me.T),1);

me.moleToMassFractions;
me.props
end