classdef SimplexTableau
    properties
        tab
        m
        n
    end
    methods
         function obj = Tableau(A, b, c)
             addpath("lib/");
             [m,n] = size(A);
             obj.m = m;
             obj.n = n;
             [B, cB, N, cN] = generate_basis(A,c,[]);
             D = B\N;
             f = B\b;
             r = cB'*D-cN';
             z = cB'*f;
             obj.tab = [zeros(1,m),r,z;eye(m),D,f];
         end
         function obj = simplex(obj)
             [~,i] = max(obj.tab(1,1:obj.n));
             d = obj.tab(2:end,i);
             if all(d <= 0)
                 error('problem unbounded. All entries <= 0 in column %d',i);
             else
                 f = obj.tab(2:end,end);
                 J = find(d>0);
                 [~,j] = min(f(J)./d(J));
                 j = J(j)+1;
                 obj.tab = gaussian_exchange(obj.tab, i, j);
                 disp(obj.tab);
             end
         end
         function [obj, optimal] = iter(obj)
             addpath("utils/");
             divider();
             if (any(obj.tab(1,1:obj.n)>0))
                 fprintf("negative r detected\n");
                 obj = obj.simplex();
                 optimal = false;
             else
                 fprintf("optimal\n");
                 optimal = true;
             end
         end
         function obj = loop(obj)
             optimal = false;
             while optimal ~= true
                 [obj, optimal] = obj.iter;
             end
         end
         function val = val(obj)
             val = obj.tab(1,end);
         end
    end
end