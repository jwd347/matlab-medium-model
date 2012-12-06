function [Teq,ATE]=equilibriumTemperature(me,Z,Tcatalyst)

ctReactants = size(me.nu',2);

% make Z be oreitned such that each column is a different mix
% and each row is a species
if size(Z,2) ~= ctReactants
    Z=Z';
end
indTest = find(abs(sum(abs(Z),2)-1)<0.0001);
ln_k_obs=me.equilibrium(Z,me.nu,me.notCondensed)';
%             indWithinT = (T>min(me.T)) & (T<max(me.T));
ctNumReactions = size(me.nu,2);

Teq=ln_k_obs.*0; % an "empty vector" of the correct dimensions
for ctReaction = 1:ctNumReactions
    Teq(indTest,ctReaction)=interp1(me.ln_kc(:,ctReaction),me.T,ln_k_obs(indTest,ctReaction));
end
if size(Tcatalyst,1)~=size(Teq,1)
    Tcatalyst=Tcatalyst';
end
ATE=repmat(Tcatalyst,1,ctNumReactions)-Teq; % ACCORDING TO MJS
if 0%any(isnan(ATE))
    figure
    hold all
    hndLines = plot(me.T-273.15,me.ln_kc,'-');
    hndSpots = plot(Tcatalyst,ln_k_obs,'.');
    disp('')
    strLines=get(hndLines);
    color = zeros(length(strLines),3);
    for i = 1:length(strLines)
        set(hndSpots(i),'Color',strLines(i).Color,'DisplayName', num2str(i,'Reaction No. %d'));
    end
    legend('Location','Best')
    title(['ATE has been found to be NaN' 10 'need to extend temperature range'])
end

end