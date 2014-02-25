function massToMoleFractions(me)

molMassInv_V = 1./me.mm_V;
me.mm = 1./(me.X*molMassInv_V);
me.Z = (me.mm*me.X)./me.mm_V';
end