%% Matthew Widjaja.
% Research Program Runner.
% Instructions: This file & func_v2.m must be in the same folder in order to be used.

clear all
clc
format short

global fixEqu
global constant
global newAlpha


%% Node Modification Query
% This lets the user select how to modify a node, if desired
fprintf('Methods: 1. AllSKO \n')
nodeMethod = input('Select a Method (1-5) = ');



%% Modification of a Node
% This lets the user change the node to a desired value

ICs = [1 1 1 1 1 1 1 1];		% Set initial conditions
fix = input('How many nodes should be fixed? (0 -> 10) = ');
	if fix >= 1
		fixEqu = zeros(1,fix);
		for i1=1:1:fix
			fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
			fixEqu(i1) = fixNode(i1);		% Saves Node # to a matrix for use in the Function
			fixValue = input('Knockout with the Value = ');		% How the node should be modified
			ICs(fixNode(i1)) = fixValue;		% Replaces the value of Node # to the declared value
		end
	else
	end
    
%% Ode45 Solver
% This solves the Diff Equs. Kinda important.

[T,Y] = ode45(@func_v2, [0 250], ICs );   

%% Create Plots
% This plots Time vs. Organism & Organism vs. Organism

graph = input('\nShould plots be generated? (y/n) = ','s');
if graph == 'y';
	yPlot = 3; xPlot = 4;		% Layout of Plots
	subplot(yPlot,xPlot,1:4);
		plot(T,Y(:,1),'-b', T,Y(:,2),'-r', T,Y(:,3),'-k', T,Y(:,4),'-g', T,Y(:,5),'--b', T,Y(:,6),'--r', T,Y(:,7),'--k', T,Y(:,8),'--g');
		legend('Org 1', 'Org 2', 'Org 3', 'Org 4', 'Org 5', 'Org 6', 'Org 7', 'Org 8');
		title('Time vs. Organism')
	subplot(yPlot,xPlot,5);
		plot(Y(:,1),Y(:,2),'-b');
		title('Org 1 vs. Org 2')
	subplot(yPlot,xPlot,6)
		plot(Y(:,2),Y(:,3),'-b')
		title('Org 2 vs. Org 3')
	subplot(yPlot,xPlot,7)
		plot(Y(:,3),Y(:,4),'-b')
		title('Org 3 vs. Org 4')
	subplot(yPlot,xPlot,8)
		plot(Y(:,4),Y(:,5),'-b')
		title('Org 4 vs. Org 5')
	subplot(yPlot,xPlot,9)
		plot(Y(:,5),Y(:,6),'-b')
		title('Org 5 vs. Org 6')
	subplot(yPlot,xPlot,10)
		plot(Y(:,6),Y(:,7),'-b')
		title('Org 6 vs. Org 1')
	subplot(yPlot,xPlot,11)
		plot(Y(:,7),Y(:,8),'-b')
		title('Org 7 vs. Org 8')
	subplot(yPlot,xPlot,12)
		plot(Y(:,8),Y(:,1),'-b')
		title('Org 8 vs. Org 1')
else
end

%% Determine Steady State Values
% This lets the user retrieve the Steady State Values

steady = input('\nShould steady-state values be printed? (y/n) = ','s');
if steady=='y'
	fprintf('The steadyState Values in the order of Node 1 --> 8 are:\n');
	for i = 1:8
		ans1(i) = Y(end,i);		% Obtains values from the last row in 'Y'
		fprintf('%f, ',ans1(i));
	end
else
end
    
%% Determine Constant & newAlpha values.
% This lets the user retrieve the Steady State Values

steady = input('\nShould constants & alpha values be printed? (y/n) = ','s');
if steady=='y'
	fprintf('The Constant Values in the order of Node 1 --> 8 are:\n');
	for i = 1:8
		ans2(i) = constant(1,i);	% Retrieves the 'Constants' matrix from func_v2.m
		fprintf('%f, ',ans2(i));
    end
    fprintf('\n\nThe Alpha Values in the order of Node 1 --> 8 are:\n');
    for i = 1:8
		ans3(i) = newAlpha(1,i);	% Retrieves the 'newAlpha' matrix from func_v2.m 
		fprintf('%f, ',ans3(i));
	end
else
end

%% Finish
% Achievement Unlocked.

fprintf('\nThis model is complete\n');