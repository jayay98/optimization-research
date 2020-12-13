n=2000;

m=200;


%%%%%%%%%%%%%%%%%%%% random data generation %%%%%%%%%%%%%%%%%%%%%
A=rand(m,n-m);

A=[A eye(m)];

b=20*rand(m,1);

c=-[rand(n-m,1);zeros(m,1)];

%%%%%%%%%%%%%%%%%%%Construct the initial simplex tableau %%%%%%%%%


tab=[-c' 0;A b];
I=n-m+1:n;

invB=eye(m,m);


%%%%%%%%%%%%%%%%%%%%%%%%%% Run simplex %%%%%%%%%%%%%%%%%%%%%%%%%


max_nb_iters=1000;
debug=0;

tic;
simplex_with_Dantzig_pivoting_rule(tab,I,m,n,max_nb_iters,debug);
toc;

%tic;
%simplex_with_Bland_pivoting_rule(tab,I,m,n,max_nb_iters,debug);
%toc;


%%%%%%%%%%%%%%%%%%%%%%%%%% Run revised simplex %%%%%%%%%%%%%%%%%%%%%%%%%

 I=n-m+1:n;
% 
 invB=eye(m,m);
% 
% 
 max_nb_iters=100;
 debug=0;
% 
 tic;
 %[c_sol,x_sol,I,invB]=revised_simplex_with_Dantzig_pivoting_rule(invB,I,m,n,A,b,c,max_nb_iters,debug);
 revised_simplex_with_Dantzig_pivoting_rule(invB,I,m,n,A,b,c,max_nb_iters,debug);
 toc;



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
 

