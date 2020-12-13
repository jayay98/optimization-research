function [matrix] = gaussian_exchange(matrix, i, j)
%i,j The pivoting element's coordinate in the matrix
    col = matrix(:, i);
    for loop = 1:length(col)
        if loop ~= j
            coeff = col(loop)/col(j);
            matrix(loop, :) = matrix(loop, :) - coeff .* matrix(j, :);
        end
    end
    matrix(j, :) = matrix(j, :) ./ col(j);
end
