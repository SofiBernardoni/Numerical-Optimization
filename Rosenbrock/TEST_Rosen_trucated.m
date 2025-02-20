%% ROSENBROCK FUNCTION
addpath("C:\Users\sofia\Documents\Numerical-Optimization")
rng(345989);

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
last_bt1=sum(btseq1)/k1; 
last_cg1=sum(cgiterseq1)/k1;

%plot
plot_iterative_optimization_results(f_Rosen,xseq1, btseq1);
figure;
hold on
plot(1:length(conv_ord1), conv_ord1, 'Color', 'b', 'LineWidth', 1.5)
title('Rosenbrock truncated Newton method [1.2, 1.2] superlinear');
hold off;

% Quadratic term of convergence
tic
[x2, f2, gradf_norm2, k2, xseq2, btseq2,cgiterseq2,conv_ord2,flag2,converged2,violations2] = truncated_newton(x0_a, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
time2=toc;
disp('Test on Rosenbrock function with x0=[1.2;1.2] and quadratic term for the adaptive tolerance: ')
disp(flag2)
disp(['Elapsed time: ',num2str(time2)])
disp(['Function value in the point found: ',num2str(f2)])
disp(['Number of violations of curvature condition: ', num2str(violations2)])
last_bt2=sum(btseq2)/k2; 
last_cg2=sum(cgiterseq2)/k2;

%plot
plot_iterative_optimization_results(f_Rosen,xseq2, btseq2);
figure;
hold on
plot(1:length(conv_ord2), conv_ord2, 'Color', 'b', 'LineWidth', 1.5)
title('Rosenbrock truncated Newton method [1.2, 1.2] quadratic');
hold off;


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
last_bt3=sum(btseq3)/k3; 
last_cg3=sum(cgiterseq3)/k3;

%plot
plot_iterative_optimization_results(f_Rosen,xseq3, btseq3);
figure;
hold on
plot(1:length(conv_ord3), conv_ord3, 'Color', 'b', 'LineWidth', 1.5)
title('Rosenbrock truncated Newton method [-1.2, 1] superlinear');
hold off;

% Quadratic term of convergence
tic
[x4, f4, gradf_norm4, k4, xseq4, btseq4,cgiterseq4,conv_ord4,flag4,converged4,violations4] = truncated_newton(x0_b, f_Rosen, gradf_Rosen, Hessf_Rosen, kmax, tolgrad, fterms_quad, cg_maxit,z0, c1, rho, btmax);
time4=toc;
disp('Test on Rosenbrock function with x0=[-1.2;1] and quadratic term for the adaptive tolerance: ')
disp(flag4)
disp(['Elapsed time: ',num2str(time4)])
disp(['Function value in the point found: ',num2str(f4)])
disp(['Number of violations of curvature condition: ', num2str(violations4)])
last_bt4=sum(btseq4)/k4; 
last_cg4=sum(cgiterseq4)/k4;

%plot
figure;
hold on
plot(1:length(conv_ord4), conv_ord4, 'Color', 'b', 'LineWidth', 1.5)
title('Rosenbrock truncated Newton method [-1.2, 1] quadratic');
hold off;
plot_iterative_optimization_results(f_Rosen,xseq4, btseq4);


%% Table with results

results_table = table({'[1.2;1.2]'; '[1.2;1.2]'; '[-1.2;1]'; '[-1.2;1]'}, ...
                      {'Superlineare'; 'Quadratica'; 'Superlineare'; 'Quadratica'}, ...
                      [f1; f2; f3; f4], [k1; k2; k3; k4], ...
                      [time1; time2; time3; time4], ...
                      [violations1; violations2; violations3; violations4],...
                      [last_cg1; last_cg2; last_cg3; last_cg4],...
                      [last_bt1; last_bt2; last_bt3; last_bt4],...
                      'VariableNames', {'Starting point',' Forcing terms', 'f_x', 'Number of iterations', 'Executing time', 'Violation of curvature conditions', 'Average of cg iterations', 'Average of bt iterations'});
writetable(results_table, 'rosenbrock_truncated.xlsx','WriteRowNames', true);
