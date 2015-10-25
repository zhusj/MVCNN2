nInstances = length(unique(imdb.images.sid));

nViews = length(imdb.images.name)/nInstances;

nDescPerShape = size(feat.x,1)/nInstances;

shapeGtClasses = imdb.images.class(1:nViews:end);

shapeSets = imdb.images.set(1:nViews:end);

nDims = size(feat.x,2);

 

% train & val

trainSets = {'train','val'};

testSets = {'test'};

[~,I] = ismember(trainSets,imdb.meta.sets);

trainSids = find(ismember(shapeSets, I));

tmp = zeros(nDescPerShape, nInstances);

tmp(:,trainSids) = 1;

trainFeat = feat.x(find(tmp)',:);

trainLabel = shapeGtClasses(trainSids)';

nTrainShapes = length(trainLabel); 

 

% test 

[~,I] = ismember(testSets,imdb.meta.sets);

testSids = find(ismember(shapeSets, I));

tmp = zeros(nDescPerShape, nInstances);

tmp(:,testSids) = 1;

testFeat = feat.x(find(tmp)',:);

testLabel = shapeGtClasses(testSids)';

nTestShapes = length(testLabel);