function a = leakyrelu(n,varargin)
%LEAKYRELU Leaky positive linear transfer function.
%	
%	See also POSLIN, PURELIN, SATLIN, SATLINS.

if nargin > 0
    n = convertStringsToChars(n);
end

if ischar(n)
  a = nnet7.transfer_fcn(mfilename,n,varargin{:});
  return
end

% Apply
a = leakyrelu.apply(n);
