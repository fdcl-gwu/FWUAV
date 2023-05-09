rng(42);
rand_idx = randperm(size(inputs, 2));
train_split = 0.8;
train_idx = rand_idx(1:round(train_split*length(rand_idx)));
val_idx = rand_idx(round(train_split*length(rand_idx)):end);
inputs_train = inputs(:, train_idx); targets_train = targets(:, train_idx);
inputs_val = inputs(:, val_idx); targets_val = targets(:, val_idx);

optimVars = [
    optimizableVariable('NHidden', [1, 3], 'Type', 'integer')
    optimizableVariable('Dropout', [0, 0.25])
%     optimizableVariable('IsDirectConnect', {'true', 'false'}, 'Type', 'categorical')
    optimizableVariable('InitialLearnRate', [3e-4, 3e-3], 'Transform', 'log')
    optimizableVariable('L2Regularization', [1e-5, 1e-3], 'Transform', 'log')
    ];

XTrain = inputs_train'; YTrain = targets_train';
XValidation = inputs_val'; YValidation = targets_val';
ObjFcn = makeObjFcn(inputs_train', targets_train', inputs_val', targets_val');

BayesObject = bayesopt(ObjFcn, optimVars, ...
    'MaxObjectiveEvaluations', 2000, ...
    'NumSeedPoints', 4, ...
    'IsObjectiveDeterministic', false, ...
    'UseParallel', true);

bestIdx = BayesObject.IndexOfMinimumTrace(end);
% fileName = BayesObject.UserDataTrace{bestIdx};
% savedStruct = load(fileName);
% ObjValue = savedStruct.ObjValue;

filename = 'iterative_learning_dnn_tuning';
save(filename, 'BayesObject');

%%
function ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation)
    ObjFcn = @PerfFun;
    N_hidden = {32, [32,96], [32,60,96]};

    function [ObjValue,cons] = PerfFun(optVars)
        drop_prob = optVars.Dropout;
        layers = [featureInputLayer(12, 'Normalization', 'rescale-symmetric')];
        for fc=N_hidden{optVars.NHidden}
            layers = [layers
                fullyConnectedLayer(fc, 'Name', "fc_"+num2str(fc))
                leakyReluLayer
                dropoutLayer(drop_prob, 'Name', "dr_"+num2str(fc))];
        end
        layers = [layers
            fullyConnectedLayer(60, 'Name', "fc_out")
            regressionLayer];
        lGraph = layerGraph(layers);

%         if optVars.IsDirectConnect == 'true'
        if true
            lGraph = disconnectLayers(lGraph, "dr_"+num2str(fc), "fc_out");
            lGraph = addLayers(lGraph, concatenationLayer(1,2,'Name','concat'));
            lGraph = connectLayers(lGraph, "dr_"+num2str(fc), 'concat/in1');
            lGraph = connectLayers(lGraph, 'input', 'concat/in2');
            lGraph = connectLayers(lGraph, 'concat/out', 'fc_out');
        end

        miniBatchSize = 6000;
        maxEpochs = 500;
        validationFrequency = 10;
        options = trainingOptions('adam', ...
            'InitialLearnRate',optVars.InitialLearnRate, ...
            'MaxEpochs',maxEpochs, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',300, ...
            'LearnRateDropFactor',0.25, ...
            'MiniBatchSize',miniBatchSize, ...
            'L2Regularization',optVars.L2Regularization, ...
            'Shuffle','every-epoch', ...
            'Verbose',false, ...
            'Plots','training-progress', ...
            'ValidationData',{XValidation,YValidation}, ...
            'ValidationFrequency',validationFrequency, ...
            'ExecutionEnvironment','cpu');

        trainedNet = trainNetwork(XTrain, YTrain, lGraph, options);
%         close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_UIFIGURE'));
        
        MSE_train = mse(trainedNet.predict(XTrain), YTrain);
        MSE_val = mse(trainedNet.predict(XValidation), YValidation);
        R2_train = regression(YTrain, trainedNet.predict(XTrain), 'one');
        R2_val = regression(YValidation, trainedNet.predict(XValidation), 'one');
        zero_err = norm(trainedNet.predict(zeros(1,12)));

%         ObjValue = zero_err + abs(R2_val - R2_train) / R2_train;
%         ObjValue = MSE_val + 3 * zero_err + 1/3 * abs(R2_val - R2_train) / R2_train;
        ObjValue = MSE_val + zero_err;
%         fileName = "bayes_hyperparameter/"+ num2str(ObjValue) + ".mat";
%         save(fileName,'trainedNet','ObjValue','options','optVars')
        cons = [];

    end
end

function [err, grad] = jac_err_weights(dlnet, inp, target, num_wts, grad_idx)
    out = predict(dlnet, inp);
    err = mse(out, target);
    grad = zeros(num_wts, 1);

    gradval = dlgradient(err, dlnet.Learnables);
    wt_idx = 0;
    for j = grad_idx
        wt_idx_new = wt_idx + numel(gradval.Value{j});
        grad((wt_idx+1):wt_idx_new) = reshape(gradval.Value{j}, [], 1);
        wt_idx = wt_idx_new;
    end
end

function [out, grad] = jac_out_weights(dlnet, inp, num_wts, grad_idx)
    out = predict(dlnet, inp);
    num_out = length(out);
    grad = zeros(num_wts, num_out);
    for i = 1:num_out
        gradval = dlgradient(out(i), dlnet.Learnables);
        wt_idx = 0;
        for j = grad_idx
            wt_idx_new = wt_idx + numel(gradval.Value{j});
            grad((wt_idx+1):wt_idx_new, i) = reshape(gradval.Value{j}, [], 1);
            wt_idx = wt_idx_new;
        end
    end
end

function weights = get_fc_weights(net, fc_idx)
    weights = [];
    for idx = fc_idx
        Layer = net.Layers(idx);
        weights = [weights; reshape(Layer.Weights, [], 1)];
        weights = [weights; Layer.Bias]; 
    end
    weights = double(weights);
end

function net = set_fc_weights(net, fc_idx, weights)
    wt_idx = 0;
    lGraph = net.layerGraph;
    for idx = fc_idx
        Layer = lGraph.Layers(idx);

        wt_idx_new = wt_idx + numel(Layer.Weights);
        Layer.Weights = reshape(weights((wt_idx+1):wt_idx_new), size(Layer.Weights));
        wt_idx = wt_idx_new;
        wt_idx_new = wt_idx + numel(Layer.Bias);
        Layer.Bias = reshape(weights((wt_idx+1):wt_idx_new), size(Layer.Bias));
        wt_idx = wt_idx_new;

        lGraph = replaceLayer(lGraph, Layer.Name, Layer);
    end

    net = assembleNetwork(lGraph);
end
