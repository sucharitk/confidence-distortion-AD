function [m,v] = beta_ab2mv(a,b)
% for a beta distribution conver mean and variance to a,b

m = a./(a+b);

v = (a.*b)./(((a+b).^2).*(a+b+1));

end