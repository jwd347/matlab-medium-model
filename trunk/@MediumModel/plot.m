function plot(me)
if length(me.T)>1   %If the plot is over a range of temperatures, produce a line graph
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
else  % if the plot is over a single temperature, produce a bar chart.

% Cp
subplot(2,1,1)
H = [me.cp_V me.cp];
N = numel(H);
names = [me.names 'Mix'];
for i=1:N
    h = bar(i, H(i));
   hold on
    col = [i/8 i/16 (-i/8)+1];
    if strcmp('H2O',names(i))==1
        col = [0 0 1]; end
      if strcmp('CO',names(i))==1
        col = [0 0 0]; end
        if strcmp('H2',names(i))==1
        col = [1 0 0]; end
        if strcmp('CO2',names(i))==1
        col = [1 1 0]; end
        if strcmp('N2',names(i))==1
            col = [0 1 0]; end
            if strcmp('O2',names(i))==1
        col = [1 1 1 ]; end    
        

    set(h, 'FaceColor', col)
end
set(gca,'Xtick',1:N,'XTickLabel',names)
ylabel('cp (J/mol K)')
xlabel('Species')

title('Graphs of stream')

% h
subplot(2,1,2)
H = [me.h_V me.h];
N = numel(H);
names = [me.names 'Mix'];
for i=1:N
    h = bar(i, H(i));
   hold on
    col = [i/8 i/16 (-i/8)+1];
     if strcmp('H2O',names(i))==1
        col = [0 0 1]; end
      if strcmp('CO',names(i))==1
        col = [0 0 0]; end
        if strcmp('H2',names(i))==1
        col = [1 0 0]; end
        if strcmp('CO2',names(i))==1
        col = [1 1 0]; end
        if strcmp('N2',names(i))==1
            col = [0 1 0]; end
            if strcmp('O2',names(i))==1
        col = [1 1 1 ]; end    
        
    set(h, 'FaceColor', col)
end
set(gca,'Xtick',1:N,'XTickLabel',names)
ylabel('h (J/mol)')
xlabel('Species')

figure
% s
subplot(2,1,1)
H = [me.s_V me.s];
N = numel(H);
names = [me.names 'Mix'];
for i=1:N
    h = bar(i, H(i));
   hold on
    col = [i/8 i/16 (-i/8)+1];
     if strcmp('H2O',names(i))==1
        col = [0 0 1]; end
      if strcmp('CO',names(i))==1
        col = [0 0 0]; end
        if strcmp('H2',names(i))==1
        col = [1 0 0]; end
        if strcmp('CO2',names(i))==1
        col = [1 1 0]; end
        if strcmp('N2',names(i))==1
            col = [0 1 0]; end
            if strcmp('O2',names(i))==1
        col = [1 1 1 ]; end    
        
    set(h, 'FaceColor', col)
end
set(gca,'Xtick',1:N,'XTickLabel',names)
ylabel('s (J/mol K)')
xlabel('Species')


%Zeq
subplot(2,1,2)
H = [me.Zeq];
N = numel(H);
names = [me.names];
for i=1:N
    h = bar(i, H(i));
   hold on
    col = [i/8 i/16 (-i/8)+1];
     if strcmp('H2O',names(i))==1
        col = [0 0 1]; end
      if strcmp('CO',names(i))==1
        col = [0 0 0]; end
        if strcmp('H2',names(i))==1
        col = [1 0 0]; end
        if strcmp('CO2',names(i))==1
        col = [1 1 0]; end
        if strcmp('N2',names(i))==1
            col = [0 1 0]; end
            if strcmp('O2',names(i))==1
        col = [1 1 1 ]; end    
        
    set(h, 'FaceColor', col)
end
set(gca,'Xtick',1:N,'XTickLabel',names)
ylabel('Composition')
xlabel('Species')


end
end