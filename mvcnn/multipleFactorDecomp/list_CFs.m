function [A,emptyCells] = list_CFs(GCF,p_datatype)
    if ~exist('p_datatype','var')
        p_datatype = 'cell';
    end
    if iscell(GCF)
        %cell array case
        A = GCF(:);
        emptyCells = cellfun(@isempty,A);
        A(emptyCells) = [];
        A = cell2mat(A')';
    elseif isstruct(GCF)
        A = {GCF.C};
        emptyCells = cellfun(@isempty,A);
        A(emptyCells) = [];
        A = cell2mat(A)';
    else %p_datatype = 'array'
        %array case
        sz = size(GCF);
        data_sz = sz(end);
        A = reshape(GCF,[],data_sz);
        
        emptyCells = zeros(size(A,1),1);
    end    
end