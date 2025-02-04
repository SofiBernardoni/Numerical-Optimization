%% ROSENBROCK FUNCTION

f_Rosen = @(x) 100*(x(2,:)-x(1,:).^2).^2+(1-x(1,:)).^2 ;
gradf_Rosen= @(x) [400*x(1,:).^3+(2-400*x(2,:)).*x(1,:)-2;
                   200*(x(2,:)-x(1,:).^2)];
Hessf_Rosen=@(x) [1200*x(1,:).^2-400*x(2,:)+2, -400*x(1,:) ;
                    -400*x(1,:), 200];
x0_a=[1.2;1.2];
x0_b=[-1.2;1];

load forcing_terms.mat

%% TEST of Trucated Newton with x0=[1.2;1.2]

kmax=500;
tolgrad=5e-7;
cg_maxit=50;

z0=zeros(2,1);
c1=1e-4;
rho=0.5;
btmax=50; % compatible with rho (with alpha0=1 you get min_step 8.8e-16)

% Superlinear term of convergence
tic
[x1, f1, gradf_norm1, k1, xseq1, btseq1,cgiterseq1,conv_ord1,flag1,converged1,violations1] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
time1=toc;
disp('Test on Rosenbrock function with x0=[1.2;1.2] and superlinear term for the adaptive tolerance: ')
disp(flag1)
disp(['Elapsed time: ',num2str(time1)])
disp(['Function value in the point found: ',num2str(f1)])
disp(['Number of violations of curvature condition: ', num2str(violations1)])
convergence_order= conv_ord1(end-5:end)
last_bt=btseq1(k1-5:k1) 
last_cg=cgiterseq1(k1-5:k1)

plot_iterative_optimization_results(f_Rosen,xseq1, btseq1);

% Quadratic term of convergence
tic
[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2,converged2,violations2] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
time2=toc;
disp('Test on Rosenbrock function with x0=[1.2;1.2] and quadratic term for the adaptive tolerance: ')
disp(flag2)
disp(['Elapsed time: ',num2str(time2)])
disp(['Function value in the point found: ',num2str(f2)])
disp(['Number of violations of curvature condition: ', num2str(violations2)])
convergence_order= conv_ord2(end-5:end)
last_bt=btseq2(k2-5:k2) 
last_cg=cgiterseq2(k2-5:k2)

plot_iterative_optimization_results(f_Rosen,xseq2, btseq2);

%% TEST of Trucated Newton  with x0=[-1.2;1]

kmax=500;
tolgrad=5e-7;
cg_maxit=50;

z0=zeros(2,1);
c1=1e-4;
rho=0.5;
btmax=50;
% rho=0.8;
% btmax=155;

% Superlinear term of convergence
tic
[x3, f3, gradf_norm3, k3, xseq3, btseq3,cgiterseq3,conv_ord3,flag3,converged3,violations3] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_suplin, cg_maxit,z0, c1, rho, btmax);
time3=toc;
disp('Test on Rosenbrock function with x0=[-1.2;1] and superlinear term for the adaptive tolerance: ')
disp(flag3)
disp(['Elapsed time: ',num2str(time3)])
disp(['Function value in the point found: ',num2str(f3)])
disp(['Number of violations of curvature condition: ', num2str(violations3)])
convergence_order= conv_ord3(end-5:end)
last_bt=btseq3(k3-5:k3) 
last_cg=cgiterseq3(k3-5:k3)

plot_iterative_optimization_results(f_Rosen,xseq3, btseq3);

% Quadratic term of convergence
tic
[x4, f4, gradf_norm4, k4, xseq4, btseq4,cgiterseq4,conv_ord4,flag4,converged4,violations4] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
time4=toc;
disp('Test on Rosenbrock function with x0=[-1.2;1] and quadratic term for the adaptive tolerance: ')
disp(flag4)
disp(['Elapsed time: ',num2str(time4)])
disp(['Function value in the point found: ',num2str(f4)])
disp(['Number of violations of curvature condition: ', num2str(violations4)])
convergence_order= conv_ord4(end-5:end)
last_bt=btseq4(k4-5:k4) 
last_cg=cgiterseq4(k4-5:k4)

plot_iterative_optimization_results(f_Rosen,xseq4, btseq4);
