%# train one-against-all models
Tt = At*MAP.R;
trainLabel = labels_train;
testLabel = test_act_labels;
numTest = length(testLabel);
%trainData = new_labels;
%testData = test_results;
trainData = MAP.T;
testData = Tt;
numLabels = 26;

trainData = gaussian_distance(MAP.T,MAP.T);
% svmStruct = svmtrain(trainLabel,trainData,'-c 26');
% C = svmclassify(svmStruct,testData,'showplot',true);
% err_rate = sum(testLabel~= C)/numTest % mis-classification rate
% conMat = confusionmat(testLabel,C)

for g = 18 %1:20
    model = ovrtrain(trainLabel, trainData, sprintf('-g %d',18));%sprintf('-c 26 -g %d -d %d -q',g));
    [pred ac decv] = ovrpredict(testLabel, testData, model);
    fprintf('Accuracy = %g%%\n', ac * 100);
    acrcy(g) = ac;
end



bestcv = 0;
for c = 26%log2c = 5 %-1:2:3,
  for g=1:20%log2g = -4:2:10,
    cmd = ['-q -c ', num2str(c), ' -g ', num2str(g)];%num2str(2^log2g)];
    cv = get_cv_ac(trainLabel, trainData, cmd, 3);
    if (cv >= bestcv),
      bestcv = cv; bestc = c; bestg =g;% 2^log2g;
    end
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', c, g, cv, bestc, bestg, bestcv);
  end
end

 model = svmtrain(double(trainLabel), trainData, '-c 26 -t 2 -b 1');
 [aa,bb,p] = svmpredict(double(testLabel), testData, model, '-b 1');

model = cell(numLabels,1);
for k=1:numLabels
    model{k} = svmtrain(double(trainLabel==k), trainData, '-c 1 -t 0 -b 1');
end


%# get probability estimates of test instances using each model
prob = zeros(numTest,numLabels);
for k=1:numLabels
    [~,~,p] = svmpredict(double(testLabel==k), testData, model{k}, '-b 1');
    prob(:,k) = p(:,model{k}.Label==1);    %# probability of class==k
end

%# predict the class with the highest probability
[~,pred] = max(prob,[],2);
acc = sum(pred == testLabel) ./ numel(testLabel)    %# accuracy
C = confusionmat(testLabel, pred) 