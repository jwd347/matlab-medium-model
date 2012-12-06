function massToMoleFractions(me)

molMassInv_V = 1./me.mm_V;
me.mm = 1./(me.X*molMassInv_V);
me.Z = (me.mm.*me.X)./molMass;
end