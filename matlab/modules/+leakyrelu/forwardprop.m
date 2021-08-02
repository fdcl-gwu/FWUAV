function da = forwardprop(dn,n,a,param)

  da = bsxfun(@times,dn,(n>=0) + 0.01*(n < 0));
end
