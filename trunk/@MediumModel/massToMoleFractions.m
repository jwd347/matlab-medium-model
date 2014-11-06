function massToMoleFractions(me)

molMassInv_V = 1./me.mm_V;
me.mm = 1./(me.X*molMassInv_V');
me.Zeq = (repmat(me.mm,1,size(me.X,2)).*me.X)./repmat(me.mm_V,size(me.X,1),1);
end