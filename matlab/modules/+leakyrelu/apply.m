function a = apply(n,param)

  a = max(0.01*n,n);
  a(isnan(n)) = nan;

end
