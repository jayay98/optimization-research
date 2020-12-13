function [B, cB, N, cN] = generate_basis(A,c,basis)
    [m,n] = size(A);
    if isempty(basis)
        basis = 1:m;
    end
    B = A(:,basis);
    cB = c(basis);
    NonBasis = find(ismember(1:n, basis)==0);
    N = A(:,NonBasis);
    cN = c(NonBasis);
end