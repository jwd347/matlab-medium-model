function moleToMassFractions(me)
me.mm = me.Zeq*me.mm_V;
me.X = me.Zeq'.*repmat(me.mm_V,1,length(me.mm))./repmat(me.mm',size(me.Zeq,2),1);
end