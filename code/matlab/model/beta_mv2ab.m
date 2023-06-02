function [a,b] = beta_mv2ab(m,v)
% for a beta distribution conver a,b to mean and variance

a = (m./v).*(m-m.*m-v);
b = (m-m.*m-v).*(1-m)./v;

end