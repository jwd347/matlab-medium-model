function setX(me,X)
if length(X)~=length(me.names)
    error('MediumModel:setX','Must have correct number of species')
end
me.X=X;
me.CheckShape;
me.massToMoleFractions;
me.props
end