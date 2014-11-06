function moleToMassFractions(me)
me.mm = me.Zeq*me.mm_V';
me.X = me.Zeq.*repmat(me.mm_V,length(me.mm),1)./repmat(me.mm,1,size(me.Zeq,2));
end