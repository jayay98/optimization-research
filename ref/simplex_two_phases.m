
%%%%%%%%%%%%%%%%% read data %%%%%%%%%%%%%%%%%%%%%

A=dlmread('config/matrix_A3.txt');
b=dlmread('config/vector_b3.txt');
c=dlmread('config/vector_c3.txt');


%A=[2 5 1 0 0; 1 3 0 -2   0; -1 -2 0 0 1;];
%b=[12;6;-5];
%c=[-1.75;-1; -0.75; -3.5;4.75];

% A=[3 1 2 4 1 0 0; 1 1 1 2 0 1 0; 1 2 3 1 0 0 1];
% b=[64;20;30];
% c=[-4; -6; -5; -1; 0;0;0];

% A=[1 1 1 2 0 0 1 0;3 1 2 4 1 -1 0 1];
% b=[20;64];
% c=-[4;6;5;3;1;-1;0;0];


% A=[1 1 1 2 0 0 1 0 0;3 1 2 4 1 -1 0 1 0;
%     1 2 3 1 1 -1 0 0 1];
%  b=[20;64;30];
%  c=-[4;6;5;3;1;-1;0;0;0];



m=size(A,1);
n=size(A,2);

%%%%%%%%%%%%%%% normalize data so that the right-hand side vector is nonnegative %%%%%%%%%%%%%%%%%

A=A.*sign(b);
b=b.*sign(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n Phase I.\n\n');
tic;
%%%%%%%%%%%%%%%%%%%%%% PHASE I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2=[A eye(m)];
c2=zeros(n+m,1);
c2(n+1:end)=ones(1,m);
I=n+1:n+m;
n2=n+m;
tab2=construct_tab(A2,b,c2,n+m,m,I);
max_nb_iters=2000;
debug=0;
[c_sol2,~,I]=simplex_with_Dantzig_pivoting_rule(tab2,I,m,n2,max_nb_iters,debug);
if(max(I)>n)
    if(abs(c_sol2)<1e-10)
        I=I(I<=n);
    else
        error('Problem is infeasible!');
    end
end

tab=construct_tab(A,b,c,n,m,I);

%%%%%%%%%%%%%%%%%%%%%% END OF PHASE I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\n Phase II.\n\n');

%%%%%%%%%%%%%%%%%%%%%% PHASE II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_nb_iters=2000;
debug=0;

[c_sol,x_sol,I,tab]=simplex_with_Dantzig_pivoting_rule(tab,I,m,n,max_nb_iters,debug);

fprintf('Simplex two Phase ');

toc;
%%%%%%%%%%%%%%%%%%%%%% END OF PHASE II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%% USE MATLAB LINEAR PROGRAM SOLVER %%%%%%%%%%%%%%%%

fprintf('\n------\n------\n \n Run Matlab solver linprog\n\n------\n------\n');
tic;
my_I=-eye(n,n);
r=zeros(n,1);
options = optimoptions('linprog','Algorithm','dual-simplex');
%options = optimoptions('linprog','Algorithm','interior-point');
lb=zeros(n,1);
ub=10^20*ones(n,1);
X=linprog(c,my_I,r,A,b,lb,ub,options);
fprintf('Linprog ');
toc;

 fprintf('The optimal objective value found by simplex is %f.\n',c'*X);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function tab=construct_tab(A,b,c,n,m,I)
if(length(I)==m)
    J=setdiff(1:n,I);
    B=A(:,I);
    N=A(:,J);
    cB=c(I);
    cN=c(J);
    if(det(B)>0)
        f=B\b;
        if(min(f)>=0)
            r=cB'*(B\N)-cN';
            t1=zeros(1,n+1);
            t1(end)=cB'*f;
            t1(J)=r;
            tab=[t1;B\A f];
        else
            error('This basic solution is not feasible!');
        end
    else
        error('B is not invertible!');
    end
else
    if(rank(A)<m)
        error('A must has rank equal to %d!', m);
    else
        J=setdiff(1:n,I);
        j=1;
        while(length(I)<m&& j<=length(J))
            if(rank(A(:,[I;j]))>rank(A(:,I)))
                I(end+1)=j;
            end
            j=j+1;
        end
        construct_tab(A,b,c,n,m,I);
    end
end

end