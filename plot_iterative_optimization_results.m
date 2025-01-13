function [flag] = plot_iterative_optimization_results(f,xseq, btseq)
% 3 plots of result of iterative method for optimization problems
% 1: top view loss f with contour and sequence generated with the method
% 2: barplot di btseq (num iterations to get steplength with backtracking)
% 3: surface view of f with surface and sequence generated with the method
% INPUT ARGUMENTS:
% f= function handle (optimization function)
% xseq = sequence of values generated with the iterative method
% xseq = sequence of num iterations to get steplength with backtracking

f_meshgrid = @(X,Y)reshape(f([X(:),Y(:)]'),size(X));

% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% Computation of the values of f for each point of the mesh
Z = f_meshgrid(X, Y);

% Plots

% Contour plot with curve levels for each point in xseq
fig1 = figure();
% ATTENTION: actually, the mesh [X, Y] is too coarse for plotting the last level curves corresponding to the last point in xseq (check it zooming
% the image).
[C2, ~] = contour(X, Y, Z, f(xseq));
hold on
% plot of the points in xseq
plot(xseq(1, :), xseq(2, :), 'r--*')
hold off

% Barplot of btseq
fig2 = figure();
bar(btseq,'r')

% Much more interesting plot
fig3 = figure();
surf(X, Y, Z,'EdgeColor','none')
hold on
plot3(xseq(1, :), xseq(2, :), f(xseq), 'r--*')
hold off

flag= 'done'; % just to have something to return (if not needed i can delete it)

end