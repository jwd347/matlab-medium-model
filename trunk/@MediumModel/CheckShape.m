function CheckShape(me)
me.T=reshape(me.T,length(me.T),1);
me.Z=reshape(me.Z,1,length(me.Z));
if size(me.nu,2)~=length(me.names) & size(me.nu,1)==length(me.names)
    me.nu = me.nu';
end
%             me.X=reshape(me.X,length(me.T),length(me.names));
%             me.T = T;
end