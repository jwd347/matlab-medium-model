function setZ(me,Z)
me.Z=Z;
me.CheckShape;
if size(me.Z,1)==1
    me.Zeq=repmat(me.Z,length(me.T),1);
else
    me.Zeq=me.Z;
end

me.props
end