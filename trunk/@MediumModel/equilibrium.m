function ln_k=equilibrium(Z,nu,notCondensed)
% ln_k has a row for each reaction.
% For many points, Z will have many rows and ln_k will have as many
% columns.

if nargin==2
    notCondensed=ones(1,size(Z,2));
end
ctReactions = size(nu,1);
ctReactants = size(nu,2);

% make Z be oriented such that each row is a different mix
% and each column is a species
if size(Z,2) ~= ctReactants
    Z=Z';
end
% for each mix cal
ctNumPoints = size(Z,1); % the number of the mix
ln_k = zeros(ctReactions,ctNumPoints);
for ctPoint = 1:ctNumPoints,
    % pick out the row of the input composition
    thisZ = Z(ctPoint,:);
    % raise each concentratrion element to the appropriate
    % power for the stoichiomtery of the reaction system
    % +ve elements of nu correspond to numerator elements
    % -ve elements of nu correspond to denominator elements
    % AND because nu has a line for each reaction
    % the Z needs copying to be the same shape.
    temp = repmat(thisZ,ctReactions,1).^(nu.*repmat(notCondensed,size(nu,1),1));
    % now take the product for each row which is the
    % equilibirum position
    k = prod(temp,2);
    
    % take the natural log
    ln_k(:,ctPoint)=log(k);
    if ~isreal(ln_k(:,ctPoint))
        ln_k(:,ctPoint)=NaN;
        %                     disp('aaa')
    end
end
end
