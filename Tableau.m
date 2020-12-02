classdef Tableau
    properties
        c
        A
        b
        basis
        
        m
        n
        
        r
        D
        f
        z
    end
    methods
        function obj = Tableau(c, A, b, basis)
            obj.c = c;
            obj.A = A;
            obj.b = b;
            obj.m = length(b);
            obj.n = length(c);
            obj.basis = basis;
        end
        function [B, cB, N, cN] = generate_basis(obj)
            B = obj.A(:,obj.basis);
            cB = obj.c(obj.basis);
            NonBasis = find(ismember(1:obj.n, obj.basis)==0);
            N = obj.A(:,NonBasis);
            cN = obj.c(NonBasis);
        end
        function tab = renew_table(obj)
            [B, cB, N, cN] = obj.generate_basis;
            obj.D = B\N;
            obj.r = cB' * obj.D - cN;
            obj.f = B\obj.b;
            obj.z = cB' * obj.f;
            tab = [zeros(1, obj.m) obj.r obj.z; eye(obj.m) obj.D obj.f];
        end
        function iterate(obj)
            
        end
    end
end