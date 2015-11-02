function new_bases = generate_new_samples(bases,dists,p_total_samples) 
debugging=0;
w_thresh = 1E-5;
len_thresh = 100;
%every row in bases is a base
%dist is a column vector by the number of bases

[dists,ordr] = sort(dists);
old_bases = bases(ordr,:);%sort the old bases
total_num_samples = 150;%this number of samples will be divided into three thirds.
if exist('p_total_samples','var')&& (p_total_samples>0)
    total_num_samples = p_total_samples;
end

%first: will be the 10% of old samples with the highest weights, or all
%samples with nonzero weights. which ever is less.
%second: samples according to weights based on old samples. the majority 80% of the total number of samples. 
%third: random combination of old samples with high weights.
%total_num_samples = length(dist);
old_num_samples = size(old_bases,1);%total_num_samples;
w = exp(-0.5*(dists.^2));
w = w - w_thresh; w(w<0)=0;
w = w/sum(w);

if nnz(w)<= 0.1*old_num_samples
    ub = nnz(w);
    %old_bases = old_bases(1:ub,:);
else
    
    % ub = ceil(length(w)*0.5);
    % while sum(w(1:ub))>0.25
    %     ub = ub*0.5;
    % end
    % ub = floor(ub*1.5);
    ub = ceil(0.10*old_num_samples);
end
new_bases = old_bases(1:ub,:);
lbls = (1:ub)';
%w = w(1:ub);
%w = w/sum(w);
% extra = total_num_samples - sum(w);
% w(1) = w(1) - extra;

num_new_samples = total_num_samples - ub;
%% generate the second part of samples
num_samples = ceil(0.90*num_new_samples);
nums = ceil(w*num_samples);
sgma = sqrt(abs(var(old_bases)));
%sgma = diag(sgma/100.0);
sgma = sgma/100.0;
%new_bases = zeros(total_num_samples ,size(bases,2));

count = 0;
for i = 1:length(w)
    if (nums(i)==0) || (count>=num_samples)% || size(new_bases,1)>old_num_samples
        break;
    end
    if count+nums(i)>num_samples
        nums(i)=num_samples-count;
    end
    %new_bases = [new_bases; mvnrnd(bases(i,:),0.25*(1-w(i))*sgma,num(i))];
    v = old_bases(i,:); 
    s = sgma;
    samples = zeros(nums(i),length(v));
    pos = 1;
    while length(v)>len_thresh 
        v1 = v(1:len_thresh);
        s1 = diag(s(1:len_thresh));
        n = mvnrnd(v1,s1/(nums(i)),nums(i));
        samples(:,pos:pos+len_thresh-1) = n;
        pos = pos +len_thresh;
        v = v(len_thresh+1:end);
        s = s(len_thresh+1:end);
    end
    n = mvnrnd(v,diag(s)/(nums(i)),nums(i));
    samples(:,pos:end) = n;
    new_bases = [new_bases; samples];
    %new_bases = [new_bases; mvnrnd(old_bases(i,:),sgma/(nums(i)),nums(i))];
    count = count + nums(i);
    
    %new_bases = [new_bases; mvnrnd(old_bases(i,:),sgma,num(i))];
    lbls = [lbls;i*ones(nums(i),1)];
end
%new_bases = new_bases(1:total_num_samples,:);%don't need this step any
%more after adding checks of count
%lbls = lbls(1:total_num_samples);

%% generate the third part of the samples.
num_samples = num_new_samples - num_samples;
facts = rand(num_samples,size(old_bases,1));
new_bases = [new_bases;facts*old_bases];

%% debug by plotting
if debugging
    fig=plotting(old_bases,(1:size(old_bases,1))');
    plotting(new_bases,lbls,1,fig,ones(length(lbls),1));
end
end