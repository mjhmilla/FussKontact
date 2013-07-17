function vec = matrix2rowVector(matrix)

rows = size(matrix,1);
cols = size(matrix,2);

vec = zeros(1, rows*cols);

for i=1:1:rows
    for j=1:1:cols
        vec( (i-1)*cols + j) = matrix(i,j);
    end    
end
