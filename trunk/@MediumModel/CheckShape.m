function CheckShape(me)
if size(me.Z,2)~=length(me.names) & size(me.Z,1)==length(me.names)
    me.Z=me.Z';
end
me.T=reshape(me.T,length(me.T),1);
if size(me.T,1)==1 & size(me.Z,1)>1
    me.T = repmat(me.T,size(me.Z,1),1);
end
if size(me.nu,2)~=length(me.names) & size(me.nu,1)==length(me.names)
    me.nu = me.nu';
end
if size(me.Z,2)~=length(me.names)
    error('MediumModel:setZ','Must have correct number of species')
end
if max(abs(sum(me.Z,2)-1)) > 1e-3
    warning('MediumModel:setZ','Composition vector does not sum to 1 (=%f)',sum(me.Z))
end

%             me.X=reshape(me.X,length(me.T),length(me.names));
%             me.T = T;
end