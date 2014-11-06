function gibbsPlot(me)
if isempty(me.ln_kc)
    error('MediumModel:gibbs','No energies calculated')
end

g_reaction = -1.*(me.R.*repmat(me.T,1,size(me.nu,1))).*me.ln_kc;
figure
plot(me.T-273.15,g_reaction,'-')
ylabel('G^o(T) [J.mol^{-1}.K^{-1}]')
xlabel('Temperature [degC]')
figure
hold all
plot(me.T-273.15,me.ln_kc,'-')
ylabel('k')
xlabel('Temperature [degC]')
end