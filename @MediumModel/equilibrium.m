function ln_k=equilibrium(Z,nu,notCondensed)
% Z must be a for a single point
if nargin==2
    notCondensed=ones(size(Z,2),1);
end
ctReactions = size(nu',1);
ctReactants = size(nu',2);

% make Z be oreitned such that each column is a different mix
% and each row is a species
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
    % AND because nu is has a line for each reaction
    % the Z needs coppyig to be the same shape.
    temp = repmat(thisZ,ctReactions,1).^(nu'.*repmat(notCondensed',size(nu,2),1));
    % now take the product for each column which is the
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
