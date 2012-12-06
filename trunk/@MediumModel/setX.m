function setX(me,X)
if length(Z)~=length(me.names)
    error('MediumModel:setZ','Must have correct number of species')
end
me.X=X;
me.CheckShape;
me.massToMoleFractions;
me.props
end