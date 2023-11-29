function param = sim_control_net(inp, IW, LW, b, inp_settings, out_settings)
%%
% param = sim(control_net, inp);

% Input process
inp = minmax_apply(inp, inp_settings{1});
inp = pca_apply(inp, inp_settings{2});
% First layer
hid = leaky_apply(IW{1} * inp + b{1});
% Second layer
inp = IW{2} * inp + LW{2,1} * hid + b{2};
% Output process
inp = minmax_reverse(inp, out_settings{1});

param = inp;

end

function a = leaky_apply(n)

  a = max(0.01*n,n);
  a(isnan(n)) = nan;

end

function y = minmax_apply(x,settings)
%MAPMINMAX.APPLY Process values

% Copyright 2012-2015 The MathWorks, Inc.

  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

function x = minmax_reverse(y,settings)
%MAPMINMAX.REVERSE Reverse process values

% Copyright 2012-2015 The MathWorks, Inc.

  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end

function y = pca_apply(x,settings)
%PROCESSPCA.APPLY Process values

  % Copyright 2012-2015 The MathWorks, Inc.

  y = settings.transform * x;
end

%% CODEGEN (not needed since direct MATLAB implementation is good enough)
% IW = control_net.IW;
% LW = control_net.LW;
% b = control_net.b;
% inp_settings = cellfun(@struct, control_net.inputs{1}.processSettings, 'UniformOutput', false);
% out_settings = cellfun(@struct, control_net.outputs{end}.processSettings, 'UniformOutput', false);

% cfg = coder.config('mex');
% cfg.InlineBetweenUserAndMathWorksFunctions = 'Always';
% codegen sim_control_net -args {inp, IW, LW, b, inp_settings, out_settings} -config cfg -launchreport;
