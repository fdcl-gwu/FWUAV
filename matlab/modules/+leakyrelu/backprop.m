function dn = backprop(da,n,a,param)

  dn = bsxfun(@times,da, (n>=0) + 0.01*(n < 0));
end
