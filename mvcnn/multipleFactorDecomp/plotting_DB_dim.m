function fig = plotting_DB_dim(DB, dim, A, L,varargin)
%plot3(A(:,1),A(:,2),A(:,3),'.','Color',[0.0700000000000001,0.680000000000000,0.0100000000000007;],'MarkerSize',1);
%plot(WL)
%colorOrder = zeros(length(A),3);
%AddMore,fig_in,estimated_L


root_folder = fullfile('C:','Amr','ChanSu','MyCode');
v = 1;
databases(v).name = 'OULUVS';
databases(v).location = fullfile(root_folder,'OULUVS');
v = v+1;
databases(v).name = 'AVLETTERS';
databases(v).location = fullfile(root_folder,'AVLETTERS');
v = v+1;
databases(v).name = '3DOBJECT';
databases(v).location = fullfile(root_folder,'3DOBJECT');

srch = strcmpi({databases.name},DB);
db_location = databases(srch).location;

load(fullfile(db_location,'configuration'));

if isnumeric(dim)
    srch = dim;
else
    srch = strcmpi({database.dimension.name},dim);
end

labels = database.dimensions(srch).values;

plotting(A, L,labels,varargin{:});

end