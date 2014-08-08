function plot(me)
if length(me.T)>1
    figure
    subplot(2,2,1), hold all
    plot(me.T,me.cp_V,'-')
    plot(me.T,me.cp,':')
    ylabel('C_p^o(T) [J.mol^{-1}.K^{-1}]')
    txLeg = {me.names{:} 'Mix' };
    legend(txLeg,'Location','Best')
    subplot(2,2,3), hold all
    plot(me.T,me.h_V,'-')
    plot(me.T,me.h,':')
    ylabel('H^o(T) [kJ.mol^{-1}]')
    txLeg = {me.names{:} 'Mix' };
    subplot(2,2,2), hold all
    plot(me.T,me.s_V,'-')
    plot(me.T,me.s,':')
    ylabel('S^o(T) [J.mol^{-1}.K^{-1}]')
    xlabel('Temperature [K]')
    subplot(2,2,4), hold all
    plot(me.T,100.*me.Zeq,'-')
    ylabel('Composition [%]')
    xlabel('Temperature [K]')
end
end