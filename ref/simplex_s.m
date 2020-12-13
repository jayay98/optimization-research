%%%%%% problem data :

 %%  A: constraint matrix  (m by n matrix)
 %%  b: right-hand side vector (n dimensional vector)
 %%  c: cost coefficient vector  (n dimensional vector)
 
 
 
%%%% tab represents the simplex tableau associated with a basic feasible
%%%% solution. tab is a m+1 by n+1 matrix.



%%%% Step 1:

if any(tab(1,1:n)>0)
    
    %%%% enter Step 2
    
    % choose pivot_column such that tab(1,pivot_column)>0
    
    d=tab(2:end,pivot_column);
    if all(d<=0)
        
        % stop and declare the unboundedness of the optimal value;
        
    else
        
        %%%% enter Step 3
        
        f=tab(2:end,end);
        
        J=find(d>0);
        
        [~,j]=min(f(J)./d(J));
        
        pivot_row=J(j)+1;
        
        %%%% enter Step 4 
        
        % apply Gauss-Jordan exchange
        
        
    end
    
    
else
    
    %%%%  stop and declare the optimality of the current solution
    
end