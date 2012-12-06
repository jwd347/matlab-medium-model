function dZ_V = reactionRate(t,Z,ln_k,nu,notCondensed)

Z=Z./sum(Z);
num_exp=abs(min(nu',0)).*repmat(notCondensed',size(nu,2),1);
den_exp = max(nu',0).*repmat(notCondensed',size(nu,2),1);
num=prod(repmat(Z',size(nu,2),1).^num_exp,2);
den = prod(repmat(Z',size(nu,2),1).^den_exp,2);
dK = num.*exp(ln_k)' - den;
%                    disp('--------------------------------')
%                    show(t,'%3.2e','t')
%                    show(dK,'%3.2e','dK')
dZ_V = (nu.*repmat(notCondensed',size(nu,2),1))*dK;
end