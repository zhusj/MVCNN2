function new_bases = resample_new_bases(bases,dist,ordr) 
%every row in bases is a base
%dist is a column vector by the number of bases
total_num_samples = 100;
bases = bases(ordr,:);
%total_num_samples = length(dist);
old_num_samples = total_num_samples;
w = exp(-0.5*(dist.^2));
w = w/sum(w);
ub = ceil(length(w)*0.5);
while sum(w(1:ub))>0.25
    ub = ub*0.5;
end
ub = floor(ub*1.5);

total_num_samples = old_num_samples - ub;
%w = w(1:ub);
%w = w/sum(w);
% extra = total_num_samples - sum(w);
% w(1) = w(1) - extra;
num = ceil(w*total_num_samples);

sgma = sqrt(abs(var(bases)));
sgma = diag(sgma/1000.0);


%new_bases = zeros(total_num_samples ,size(bases,2));
new_bases = bases(1:ub,:);
for i = 1:length(w)
    if num(i)==0 || size(new_bases,1)>old_num_samples
        break;
    end
    %new_bases = [new_bases; mvnrnd(bases(i,:),0.25*(1-w(i))*sgma,num(i))];
    new_bases = [new_bases; mvnrnd(bases(i,:),sgma.^(num(i)),num(i))];
end
    
new_bases = new_bases(1:old_num_samples,:);

