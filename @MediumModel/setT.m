function setT(me,T)
me.T=T;
me.CheckShape;
me.setZ(me.Z); % sets me.Zeq and calls moleToMassFractions & props
end