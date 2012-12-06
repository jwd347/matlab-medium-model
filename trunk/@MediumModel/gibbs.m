function g_reaction=gibbs(me)
% GIBBS Calculate the Gibbs Free energy
% for the reaction system specified in preoperty "nu"
if isempty(me.nu)
    error('MediumModel:gibbs','No reaction system defined')
end
g_V=me.h_V-(me.s_V.*repmat(me.T,1,length(me.names))); %#ok<PROP>
g_reaction=zeros(size(me.T,1),size(me.nu,2));
for ctT = 1:length(me.T)
    g_reaction(ctT,:) = g_V(ctT,:)*(me.nu); % compute standard enthalpy change of reaction
end

log_k=-g_reaction./(me.R.*repmat(me.T,1,size(me.nu,2)));
me.ln_kc=log_k;
end