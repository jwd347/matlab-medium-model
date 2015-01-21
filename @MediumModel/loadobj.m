function obj = loadobj(me)
% Function to transpose properties if they were saved the wrong way around.
% These properties were changed in an earlier update.

if size(me.Z,2)~=length(me.names)
    me.mm_V = me.mm_V';
    me.nu = me.nu';
    me.Z = me.Z';
    me.X = me.X';
    me.notCondensed = me.notCondensed';
end
obj = me;
