function write_matrix_for_latex(A)
    diary 'matlab_resutls.txt';
    fprintf('\\begin{bmatrix}\n');
    for i = 1:size(A,1)
        for j = 1:size(A,2)-1
            fprintf('%.1f & ',A(i,j));
        end
        fprintf('%.1f\\\\ \n',A(i,end));
    end
    fprintf('\\end{bmatrix}\n');
    diary off
end