function g_reaction=gibbs(me)
% GIBBS Calculate the Gibbs Free energy
% for the reaction system specified in preoperty "nu"
if isempty(me.nu)
    error('MediumModel:gibbs','No reaction system defined')
end
g_V=me.h_V-(me.s_V.*repmat(me.T,1,length(me.names))); %#ok<PROP>
g_reaction=zeros(size(me.T,1),size(me.nu,2));
g_reaction = g_V*me.nu'; % compute standard enthalpy change of reaction

log_k=-g_reaction./(me.R.*repmat(me.T,1,size(me.nu,1)));
me.ln_kc=log_k;
end