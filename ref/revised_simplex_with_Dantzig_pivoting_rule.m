function [c_sol,x_sol,I,invB]=revised_simplex_with_Dantzig_pivoting_rule(invB,I,m,n,A,b,c,max_nb_iters,debug)

validate_input(invB,I,m,A,b); % check if B is a basis matrix associated with a basic feasible solution

tmp=1:n;

keep_running=true;
nb_iters=0;

while keep_running&& nb_iters<max_nb_iters
    nb_iters=nb_iters+1;
    T=setdiff(tmp,I);  % get indices of nonbasic variables
    r=c(I)'*invB*A(:,T)-c(T)';  % compute the reduced cost coefficient
    if any(r>0) %check if there is positive reduced cost coefficient
        [~,tmp_k] = max(r); %yes, find the largest one (Dantzig's pivoting rule)
        pivot_column=T(tmp_k);
        d=invB*A(:,pivot_column);  % compute the column vector
        if all(d<=0) %check if the corresponding column is negative
            error('problem unbounded. All entries <= 0 in column %d',pivot_column); % yes, the problem is unbounded
        else %otherwise, find pivoting row  
            f=invB*b;
            J=find(d>0); %find the indices corresponding to positive entries of d;
            [~,j]=min(f(J)./d(J));
            pivot_row=J(j);
            if debug==1
                fprintf('pivot row is %d\n',I(pivot_row));
                fprintf('pivot column is %d\n',pivot_column);
            end
   
            %Apply Gauss-Jordan Exchange
            
            D=[invB d];
            
            %normalize
            D(pivot_row,:) = D(pivot_row,:)/D(pivot_row,end);    
            %construct the Q matrix
            d_full=-D(:,end);
            d_full(pivot_row)=1;
            Q=eye(m,m);
            Q(:,pivot_row)=d_full;
            %do elementary row operations
            D=Q*D;
            %update the basis index
            I(pivot_row)=pivot_column;
            %update the inverse of B
            invB=D(1:m,1:m);
        end
        if debug==1  %print current basic feasible solution
            fprintf('Simplex Iteration %d \n : current basic feasible solution is\n',nb_iters);
            disp(get_current_x(invB,I,b,n));
        end
      %  fprintf('Simplex Iteration %d \n : current objective value is %f.\n',nb_iters,c'*get_current_x(invB,I,b,n));
    else
        keep_running=false;
        fprintf('Optimal Solution found in %d simplex iterations.\n',nb_iters);
        if debug==1  %print current basic feasible solution
            fprintf('The optimal basic feasible solution is:\n');
            disp(get_current_x(invB,I,b,n));
        end
        fprintf('The optimal objective value found by simplex is %f.\n',c'*get_current_x(invB,I,b,n));
    end
end
x_sol=get_current_x(invB,I,b,n);
c_sol=c'*x_sol;
end


function current_x = get_current_x(invB,I, b,n)
current_x = zeros(n,1);
current_x(I)=invB*b;
end

function validate_input(invB,I,m,A,b)
if ~ismatrix(invB)
    error('B must be a matrix');
end
if size(invB,1)~=size(invB,2)
    error('invB must be a square matrix');
end
if ~(max(max(abs(A(:,I)*invB-eye(m))))<1e-10)
    error('invB must be the inverse of B');
end
if min(invB*b)<0
    error('invB does not give a feasible solution');
end
end