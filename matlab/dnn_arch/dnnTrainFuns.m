function funs = dnnTrainFuns
  funs.norm_data=@norm_data;
  funs.train_val_ds=@train_val_ds;
  funs.initialize_nn_param=@initialize_nn_param;
  funs.initializeGlorot=@initializeGlorot;
  funs.train_custom_network=@train_custom_network;
  funs.modelLoss=@modelLoss;
  funs.forwardnet=@forwardnet;
  funs.predictnet=@predictnet;
  funs.thresholdL2Norm=@thresholdL2Norm;

  funs.get_dl_weights = @get_dl_weights;
  funs.set_dl_weights = @set_dl_weights;
  funs.jac_err_weights = @jac_err_weights;
  funs.jac_out_weights = @jac_out_weights;
end

function [nn_inp1, nn_inp2, nn_out] = norm_data(data_inp1, data_inp2, data_out, ...
    inputs_mean, inputs_scale, targets_mean, targets_scale)
    nn_inp1 = (dlarray(data_inp1, "CB") - inputs_mean{1}) ./ inputs_scale{1};
    nn_inp2 = (dlarray(data_inp2, "CB") - inputs_mean{2}) ./ inputs_scale{2};
    nn_out = (dlarray(data_out, "CB") - targets_mean) ./ targets_scale;
end

function [mbq, mbq_val, mbq_full] = train_val_ds(inputs_multiple, targets, split_ratio, miniBatchSize, OutputEnvironment)
    N_data = size(targets, 2);
    rand_idx = randperm(N_data);
    train_idx = rand_idx(1:round(split_ratio*N_data));
    val_idx = rand_idx((round(split_ratio*N_data)+1):end);

    dsTrain = combine(arrayDatastore(inputs_multiple{1}(:, train_idx)'), ...
        arrayDatastore(inputs_multiple{2}(:, train_idx)'), ...
        arrayDatastore(targets(:, train_idx)'));
    dsVal = combine(arrayDatastore(inputs_multiple{1}(:, val_idx)'), ...
        arrayDatastore(inputs_multiple{2}(:, val_idx)'), ...
        arrayDatastore(targets(:, val_idx)'));

    mbq = minibatchqueue(dsTrain,...
        MiniBatchSize=miniBatchSize,...
        MiniBatchFormat={'BC','BC','BC'},...
        PartialMiniBatch="discard", ...
        OutputEnvironment=OutputEnvironment);
    mbq_full = minibatchqueue(dsTrain,...
        MiniBatchSize=length(train_idx),...
        MiniBatchFormat={'BC','BC','BC'},...
        PartialMiniBatch="discard", ...
        OutputEnvironment=OutputEnvironment);
    if ~isempty(val_idx)
        mbq_val = minibatchqueue(dsVal,...
            MiniBatchSize=length(val_idx),...
            MiniBatchFormat={'BC','BC','BC'},...
            PartialMiniBatch="discard", ...
            OutputEnvironment=OutputEnvironment);
    else
        mbq_val = [];
    end
end

function [dlnet, netParams] = initialize_nn_param(layers)
    dlnet = dlnetwork(layers);
    numHidFeat = numel(dlnet.Learnables.Value{end});

%     idx = 1 : length(dlnet.Learnables.Value);
    idx = dlnet.Learnables.Parameter == "Bias";
    dlnet.Learnables(idx, :) = dlupdate(@(wb) 0.01 * initializeRands(wb.size), dlnet.Learnables(idx, :));

    netParams = {
%         "conv_1d","Weights", {dlarray(0.1*sqrt(2/(1+2))*randn(1, 2, 'single'))};
%         "conv_1d","Bias", {dlarray(0.00*randn(1, 1, 'single'))};
        "conv_1d","Weights", {initializeGlorot([numHidFeat, 2], numHidFeat, numHidFeat)};
%         "conv_1d","Bias", {dlarray(zeros([numHidFeat, 1], 'single'))};
        "conv_1d","Bias", {0.01 * initializeRands([numHidFeat, 1])};
        "fc_out","Weights", {dlarray(initializeGlorot([60, numHidFeat+12], 60, numHidFeat+12))};
%         "fc_out","Bias", {dlarray(zeros([60, 1], 'single'))}
        "fc_out","Bias", {0.01 * initializeRands([60, 1])}
        };
    netParams = cell2table(netParams, 'VariableNames', {'Layer', 'Parameter', 'Value'});

end

function weights = initializeRands(sz)
weights = single(rands(prod(sz)));
weights = dlarray(reshape(weights, sz));

end

function weights = initializeGlorot(sz,numOut,numIn,className)
arguments
    sz
    numOut
    numIn
    className = 'single'
end

Z = 2*rand(sz,className) - 1;
bound = sqrt(6 / (numIn + numOut));
weights = bound * Z;
weights = dlarray(weights);

end

function [dlnet, netParams, loss] = train_custom_network(dlnet, netParams, mbq, ...
    numEpochs, learningRate, executionEnvironment, ...
    mbq_val, validationFrequency, validationCheck, ...
    regularizationCoeff)

% fig = figure;
% lineLossTrain = animatedline(Color=[0 0.4470 0.7410]);
% lineLossVal = animatedline(Color=[0.8500 0.3250 0.0980]);
% ylim([0 inf])
% xlabel("Iteration")
% ylabel("Loss")
% grid on

% fig = figure;
% ax_w = subplot(2,1,1);
% ax_g = subplot(2,1,2);

monitor = trainingProgressMonitor;
monitor.Metrics = ["TrainingLoss","ValidationLoss","Conv1DBias","Conv1DWeight1","Conv1DWeight2"];
groupSubPlot(monitor,"Loss",["TrainingLoss","ValidationLoss"]);
% groupSubPlot(monitor,"R2",["ValidationR2"]);
groupSubPlot(monitor,"Conv1D",["Conv1DWeight1","Conv1DWeight2","Conv1DBias"]);
monitor.Info = "Epoch";
monitor.XLabel = "Iteration";
monitor.Progress = 0;

trailingAvgSubnet = []; trailingAvgSqSubnet = [];
trailingAvgParams = []; trailingAvgSqParams = [];
gradDecay = 0.9;
gradDecaySq = 0.999;
regularizationDecay = exp(log(1e-4) / numEpochs); % Number is final/initial

start = tic;
iteration = 0;
% Loop over mini-batches.
% for epoch = 1:numEpochs
epoch = 0;
while epoch < numEpochs && ~monitor.Stop
    epoch = epoch + 1;

%     regularizationCoeff = regularizationCoeff_init * exp((epoch-1) * log(1e-2) / numEpochs); % exponential decay
%     regularizationCoeff = regularizationCoeff * regularizationDecay;
%     if mod(epoch, round(0.5 * numEpochs)) == 0
%         learningRate = 0.1 * learningRate;
%     end

    % Reset and shuffle mini-batch queue
    shuffle(mbq);

    while hasdata(mbq)
        iteration = iteration + 1;
        [inp1, inp2, out] = next(mbq);

        % If training on a GPU, then convert data to gpuArray.
        if executionEnvironment == "gpu"
            inp1 = gpuArray(inp1); inp2 = gpuArray(inp2); out = gpuArray(out);
        end
    
        % Evaluate the model loss and gradients using dlfeval and the modelLoss
        [loss,gradientsSubnet,gradientsParams] = dlfeval(@modelLoss,dlnet,netParams,inp1,inp2,out,regularizationCoeff);
        % L2 Regularization
        gradientThreshold = 1e1;
        idx = dlnet.Learnables.Parameter == "Weights";
%         gradientsSubnet(idx,:) = dlupdate(@(g,w) g + regularizationCoeff*w, gradientsSubnet(idx,:), dlnet.Learnables(idx,:));
        gradientsSubnet = dlupdate(@(g) thresholdL2Norm(g, gradientThreshold), gradientsSubnet);
        idx = (netParams.Parameter == "Weights");% & (netParams.Layer ~= "conv_1d");
%         gradientsParams(idx,:) = dlupdate(@(g,w) g + regularizationCoeff*w, gradientsParams(idx,:), netParams(idx,:));
        gradientsParams = dlupdate(@(g) thresholdL2Norm(g, gradientThreshold), gradientsParams);
    
        % Update the Siamese subnetwork parameters.
        [dlnet,trailingAvgSubnet,trailingAvgSqSubnet] = adamupdate(dlnet,gradientsSubnet, ...
            trailingAvgSubnet,trailingAvgSqSubnet,iteration,learningRate,gradDecay,gradDecaySq);
%         [dlnet,trailingAvgSubnet] = sgdmupdate(dlnet,gradientsSubnet, ...
%             trailingAvgSubnet,learningRate,gradDecay);
    
        % Update the fullyconnect parameters.
        [netParams,trailingAvgParams,trailingAvgSqParams] = adamupdate(netParams,gradientsParams, ...
            trailingAvgParams,trailingAvgSqParams,iteration,learningRate,gradDecay,gradDecaySq);
%         [netParams,trailingAvgParams] = sgdmupdate(netParams,gradientsParams, ...
%             trailingAvgParams,learningRate,gradDecay);
    
        % Update the training loss progress plot.
%         D = duration(0,0,toc(start),Format="hh:mm:ss");
%         lossValue = double(loss);
%         addpoints(lineLossTrain,iteration,lossValue);
%         title("Elapsed: " + string(D))
%         drawnow limitrate
        recordMetrics(monitor,iteration, ...
            TrainingLoss=loss);
        updateInfo(monitor,Epoch=string(epoch) + " of " + string(numEpochs));
        monitor.Progress = 100*epoch/numEpochs;
    end

%     if mod(epoch, 5) == 0
%         gradnet.Learnables = gradientsSubnet;
%         gradLoss = get_dl_weights(gradnet, gradientsParams, []);
%         weights = get_dl_weights(dlnet, netParams, []);
%         plot(ax_w, weights)
%         title(ax_w, "Epoch " + string(epoch))
%         plot(ax_g, gradLoss)
%         title(ax_g, "Regularization " + string(regularizationCoeff) + " LearningRate " + string(learningRate))
%     end

    if (mod(epoch, validationFrequency) == 0) && (validationCheck)
        shuffle(mbq_val);
        [inp1_val, inp2_val, out_val] = next(mbq_val);
        Y_val = predictnet(dlnet, netParams, inp1_val, inp2_val);
        ValLoss = double(l2loss(Y_val, out_val));
%         addpoints(lineLossVal,iteration,ValLoss);
%         drawnow limitrate
%         ValR2 = regression(extractdata(out_val)', extractdata(Y_val)', 'one');
        recordMetrics(monitor,iteration, ...
            ValidationLoss=ValLoss, ...
            Conv1DWeight1=netParams.Value{1}(1), ...
            Conv1DWeight2=netParams.Value{1}(2), ...
            Conv1DBias=netParams.Value{2}(1));
    end
end

% close(fig);
end

function [loss,gradientsSubnet,gradientsParams] = modelLoss(dlnet, netParams, inp1, inp2, out, regularizationCoeff)
    Y = forwardnet(dlnet, netParams, inp1, inp2);
%     loss = mse(Y, out);

    l2regularization = 0;
    idx = 1:length(dlnet.Learnables.Value);
    idx = idx(dlnet.Learnables.Parameter == "Weights");
    for i=idx
        % Spectral norm regularization
        l2regularization = l2regularization + norm(dlnet.Learnables.Value{i}.extractdata)^2;
%         sum(dlnet.Learnables.Value{i}.^2, 'all');
    end
    idx = 1:length(netParams.Value);
    idx = idx(netParams.Parameter == "Weights");
    for i=idx
        l2regularization = l2regularization + norm(netParams.Value{i}.extractdata)^2; 
    end

    loss = l2loss(Y, out) + regularizationCoeff* l2regularization;
    [gradientsSubnet,gradientsParams] = dlgradient(loss,dlnet.Learnables,netParams);
end

function Y = forwardnet(dlnet, netParams, inp1, inp2)
    Y1 = forward(dlnet, inp1);
    Y2 = forward(dlnet, inp2);

%     Y = cat(3, Y1, Y2);
%     Y = dlarray(permute(stripdims(Y), [3,1,2]), 'CBT');
%     Y = Y1 * netParams.Value{1}(1) + Y2 * netParams.Value{1}(2) + netParams.Value{2};
    Y = sum(cat(3, Y1, Y2) .* dlarray(netParams.Value{1}, "CU"), 3) + netParams.Value{2};
%     Y = leakyrelu(Y);

    Y = cat(1, Y, inp2);
    Y = fullyconnect(Y,netParams.Value{3},netParams.Value{4});
end

function Y = predictnet(dlnet, netParams, inp1, inp2)
    Y1 = predict(dlnet, inp1);
    Y2 = predict(dlnet, inp2);

%     Y = cat(3, Y1, Y2);
%     Y = dlarray(permute(stripdims(Y), [3,1,2]), 'CBT');
%     Y = Y1 * netParams.Value{1}(1) + Y2 * netParams.Value{1}(2) + netParams.Value{2};
    Y = sum(cat(3, Y1, Y2) .* dlarray(netParams.Value{1}, "CU"), 3) + netParams.Value{2};
%     Y = leakyrelu(Y);

    Y = cat(1, Y, inp2);
    Y = fullyconnect(Y,netParams.Value{3},netParams.Value{4});
end

function gradients = thresholdL2Norm(gradients,gradientThreshold)
% Gradient Clipping
gradientNorm = sqrt(sum(gradients(:).^2));
if gradientNorm > gradientThreshold
    gradients = gradients * (gradientThreshold / gradientNorm);
end

end

function weights = get_dl_weights(dlnet, netParams, fc_idx)
    weights = [];
    for idx = 1:length(dlnet.Learnables.Value)
        weights = [weights; reshape(dlnet.Learnables.Value{idx}, [], 1)];
    end
    for idx = 1:length(netParams.Value)
        weights = [weights; reshape(netParams.Value{idx}, [], 1)];
    end
    weights = double(extractdata(weights));
end

function [dlnet, netParams] = set_dl_weights(dlnet, netParams, fc_idx, weights)
    wt_idx = 0;
    weights = dlarray(weights);
    for idx = 1:length(dlnet.Learnables.Value)
        Value = dlnet.Learnables.Value{idx};
        wt_idx_new = wt_idx + numel(Value);
        dlnet.Learnables.Value{idx} = single(reshape(weights((wt_idx+1):wt_idx_new), size(Value)));
        wt_idx = wt_idx_new;
    end
    for idx = 1:length(netParams.Value)
        Value = netParams.Value{idx};
        wt_idx_new = wt_idx + numel(Value);
        netParams.Value{idx} = single(reshape(weights((wt_idx+1):wt_idx_new), size(Value)));
        wt_idx = wt_idx_new;
    end
end

function [err, grad] = jac_err_weights(dlnet, netParams, inp1, inp2, target, num_wts, grad_idx)
    out = predictnet(dlnet, netParams, inp1, inp2);
    err = mse(out, target);
    grad = zeros(num_wts, 1);

    [gradSubnet, gradParams] = dlgradient(err, dlnet.Learnables, netParams);
    wt_idx = 0;
    for gradval={gradSubnet, gradParams}
        Value = gradval{1}.Value;
        for j = 1:length(Value)
            wt_idx_new = wt_idx + numel(Value{j});
            grad((wt_idx+1):wt_idx_new) = reshape(Value{j}, [], 1);
            wt_idx = wt_idx_new;
        end
    end
end

function [out, grad] = jac_out_weights(dlnet, netParams, inp1, inp2, num_wts, grad_idx)
    out = predictnet(dlnet, netParams, inp1, inp2);
    num_out = length(out);
    grad = zeros(num_wts, num_out);
    for i = 1:num_out
        [gradSubnet, gradParams] = dlgradient(out(i), dlnet.Learnables, netParams);
        wt_idx = 0;
        for gradval={gradSubnet, gradParams}
            Value = gradval{1}.Value;
            for j = 1:length(Value)
                wt_idx_new = wt_idx + numel(Value{j});
                grad((wt_idx+1):wt_idx_new, i) = reshape(Value{j}, [], 1);
                wt_idx = wt_idx_new;
            end
        end
    end
end
