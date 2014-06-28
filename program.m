%% Matthew Widjaja.
% Research Program Runner.

% Instructions: This & func.m must be in the same folder &
% this folder must be the 'Current Folder' in MATLAB.

clear all
clc
format short
global fixEqu

%% Modification of a Node
% This lets the user change the node to a desired value

ICs = [1 1 1 1 1 1 1 1 1 1];   % Set initial conditions
fix = input('How many nodes should be fixed? (0 -> 10) = ');
	if fix >= 1
		fixEqu = zeros(1,fix);
		for i=1:1:fix
			fixNode(i) = input('Knockout Node Number = ');   % Node # to modify
			fixEqu(i) = fixNode(i);   % Saves Node # to a matrix for use in the Function
			fixValue = input('Knockout with the Value = ');   % How the node should be modified
			ICs(fixNode(i)) = fixValue;   % Replaces the value of Node # to the declared value
		end
	else
	end
    
%% Ode45 Solver
% This solves the Diff Equs. Kinda important.

[T,Y] = ode45(@func, [0 1000], ICs );   

%% Create Plots
% This plots Time vs. Organism & Organism vs. Organism

graph = input('\nShould plots be generated? (y/n) = ','s');
if graph == 'y';
	xPlot = 4; yPlot = 5;  % Layout of Plots
	subplot(xPlot,yPlot,1:10);
		plot(T,Y(:,1),'-b', T,Y(:,2),'-r', T,Y(:,3),'-k', T,Y(:,4),'-g', T,Y(:,5),'-c', T,Y(:,6),'--b', T,Y(:,7),'--r', T,Y(:,8),'--k', T,Y(:,9),'--g', T,Y(:,10),'--c');
		legend('Org 1', 'Org 2', 'Org 3', 'Org 4', 'Org 5', 'Org 6', 'Org 7', 'Org 8', 'Org 9', 'Org 10');
		title('Time vs. Organism')
	subplot(xPlot,yPlot,11);
		plot(Y(:,1),Y(:,2),'-b');
		title('Org 1 vs. Org 2')
	subplot(xPlot,yPlot,12)
		plot(Y(:,2),Y(:,3),'-b')
		title('Org 2 vs. Org 3')
	subplot(xPlot,yPlot,13)
		plot(Y(:,3),Y(:,4),'-b')
		title('Org 3 vs. Org 4')
	subplot(xPlot,yPlot,14)
		plot(Y(:,4),Y(:,5),'-b')
		title('Org 4 vs. Org 5')
	subplot(xPlot,yPlot,15)
		plot(Y(:,5),Y(:,6),'-b')
		title('Org 5 vs. Org 6')
	subplot(xPlot,yPlot,16)
		plot(Y(:,6),Y(:,7),'-b')
		title('Org 6 vs. Org 7')
	subplot(xPlot,yPlot,17)
		plot(Y(:,7),Y(:,8),'-b')
		title('Org 7 vs. Org 8')
	subplot(xPlot,yPlot,18)
		plot(Y(:,8),Y(:,9),'-b')
		title('Org 8 vs. Org 9')
	subplot(xPlot,yPlot,19)
		plot(Y(:,9),Y(:,10),'-b')
		title('Org 9 vs. Org 10')
	subplot(xPlot,yPlot,20)
		plot(Y(:,10),Y(:,1),'-b')
		title('Org 10 vs. Org 1')
else
end

%% Save Matrix Query
% This lets the user save the T & Y Matrix to a Spreadsheet

file = input('\nShould the T vs. Y Matrix be saved to a spreadsheet? (y/n) = ','s');
if file=='y'
	fileName = input('Enter in the File Name w/o the .csv extension = ','s');
	fileExt = '.csv';
	fileName = strcat(fileName, fileExt);   % Adds .csv to the desired fileName
	fileMatrix = [T Y];   % Creates a Matrix with the values of T & Y
	csvwrite(fileName,fileMatrix);
else
end

%% Determine Steady State Values
% This lets the user determine the Steady State Values

steady = input('\nShould steady-state values be calculated? (y/n) = ','s');
if steady=='y'
	i = 11;   % This confirms that the while loop does not end early
	steadyNode = input('Which node should determine Steady State? (#) = ');
	while i <= numel(T);   % While 'i' is less than the Max-Time interval tested
		if abs(Y(i,steadyNode) - Y(i-10,steadyNode)) <= 0.00005   % This tests for Steady State
			fprintf('Node %g reached a steady state\n', steadyNode);
			for j=1:10
				steadyValue(j,1) = j;   % This makes column 1 be values of 1 to 10
				steadyValue(j,2) = Y(1000,j);   % This makes column 2 be the equilibrium reached
				fprintf('Node %g \t %g \n',steadyValue(j,1),steadyValue(j,2));
			end
			i = numel(T) + 1;
		else
			i = i + 1;
		end
	end
else
end
    
%% Finish
% Achievement Unlocked.

fprintf('\nThis model is complete\n');