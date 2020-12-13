%%%%%% problem data :

 %%  A: constraint matrix  (m by n matrix)
 %%  b: right-hand side vector (n dimensional vector)
 %%  c: cost coefficient vector  (n dimensional vector)
 
 
 
%%%% invB represents the inverse of the current basis matrix B. invB
%%is an m by m matrix.



%%%% Step 1:

% compute r

if any(r>0)
    
    %%%% enter Step 2
    
    % choose pivot_column such that tab(1,pivot_column)>0
    
    d=invB*A(:,pivot_column);
    
    if all(d<=0)
        
        % stop and declare the unboundedness of the optimal value;
        
    else
        
        %%%% enter Step 3
        
        f=invB*b;
        
        J=find(d>0);
        
        [~,j]=min(f(J)./d(J));
        
        pivot_row=J(j);
        
        %%%% enter Step 4 and 5 
        
        % compute the new invB
        
        
    end
    
    
else
    
    %%%%  stop and declare the optimality of the current solution
    
end