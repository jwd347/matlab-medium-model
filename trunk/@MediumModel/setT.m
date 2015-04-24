function setT(me,T)
me.T=T;
me.CheckShape;
me.setZ(me.Z); % sets me.Zeq to match size(T) and calls moleToMassFractions & props
end