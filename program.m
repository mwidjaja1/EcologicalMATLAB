%% Matthew Widjaja.
% Research Program Runner.
% Instructions: This file & func_v4.m must be in the same folder in order to be used.

clear all
format short

global equName
global alteredIC
global fixEqu
global maxNode
global newAlpha
global steadyValueP
global oldAlpha
global constantR


%% Blank Slate
% Query to clear Command Window
userQuery = input('Should Command Window be cleared? (y/n) = ','s');
if userQuery ~= 'n'
	clc
end


%% General Parameters
% General Parameters to help model efficiency
maxNode = 8;		% Amount of Nodes present
masterNode = 2;		% Amount of Master Nodes present
	masterNodes = [4 8];	% States which nodes are Masters -- For Method 8
maxTime = 200;		% Max amount of time to use
equName = @func_v4;		% Name of file w. Equations
masterIC = ones(1,maxNode);		% Set master initial conditions
alteredIC = masterIC;		% Creates a matrix for the altered ICs
intCount = 1;		% Integer Counter for some methods


%% General Data
% General Data to run the model
steadyValueP = [38.8089080995067, 43.5963099130963, 39.4611963199623, ...
 27.6303237065429, 26.8372058084276, 28.0030748854752, ...
 34.7552952146415, 47.0617038749445];	
oldAlpha = [-2.3305, -4.7368, -0.8303, -1.4475, 0, 0, 0, -0.291868997002398;  -0.2565, -4.7009, -0.522, -0.3013, 0, 0, 0, -0.647544780938589; -0.1342, -1.1164, -1.7607, -0.5695, 0, 0, 0, -0.523688100581974; -0.2881, -0.3794, -0.2665, -2.968, 0, 0, 0, -1.01620423202169;  0, 0, 0, -0.478736026755033, -3.1091, -0.0561, -1.0181, -0.8356; 0, 0, 0, -1.08243797034983, -0.0955, -2.1481, -1.1953, -0.347; 0, 0, 0, -0.242649494404053, -0.8154, -0.4252, -2.7833, -0.576; 0, 0, 0, -1.12489771091476, -0.3553, -0.0653, -0.3707, -2.3896];						


%% Method 8 Specific Parameters
% These parameters set different weights & goals for each node in Method 8
weightValue = 0.5;
for i1 = 1:maxNode
	weightValues(i1) = weightValue;
end
% weightValues = [0, 0, 1, 0, 1, 0, 0, 0];	
desiredValues = [18.9308, 17.9217, 0, 57.0000, 45.7847, 12.9408, 0.0010, 10.0000]; 


%% Node Modification Query
% This lets the user select how to modify a node, if desired
fprintf('\nPossible Methods:');
fprintf('\n1. SKO for all Nodes');
fprintf('\n2. DKO for all Nodes');
fprintf('\n3. TKO for Nodes 4 & 8');
fprintf('\n4. Scaling Master Nodes'); 
fprintf('\n5. Fix Nodes');
fprintf('\n6. Sweep 1 Node');
fprintf('\n7. Sweep 2 Nodes');
fprintf('\n8. Manipulate Slave Nodes for a Master Node');

fprintf('\n\n');
nodeMethod = input('Select a Method (0-8) = ');


%% All Methods -- WT Calculation
% This calculates the WT Values for a second Results Matrix called WT
intCount = intCount - 1;
fixEqu = 0;
[T,W] = ode45(equName, [0 maxTime], masterIC );		% Solves the model
WT(1,:) = W(end,:);		% Saves WT data


%% Method 1: SKO for All -- Solver
% Solves when all nodes get single knocked out
if nodeMethod == 1
	for i1 = 1:maxNode
		fixEqu = i1;	% States equation to fix
		alteredIC(i1) = 0;	% Fixes the node's IC to 0
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		yData(intCount,:) = Y(end,:);	% Saves data
		xData(intCount,:) = i1;
		alteredIC = masterIC;	% Resets the IC to default values
		intCount = intCount + 1;
	end
end
	
	
%% Method 2: DKO for All -- Solver
% Solves when all nodes get double knocked out
if nodeMethod == 2
	for i1 = 1:maxNode-1
		fixEqu(1) = i1;		% States equation to fix
		alteredIC(i1) = 0;	% Fixes the node's IC to 0
		for i2 = (i1+1):maxNode
			fixEqu(2) = i2;		% States equation to fix
			alteredIC(i2) = 0;	% Fixes the node's IC to 0
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			yData(intCount,:) = Y(end,:);	% Saves data
			xData(intCount,:) = [i1, i2];
			alteredIC(i2) = masterIC(i2);	% Resets the second looping IC to default
			intCount = intCount + 1;
		end
		alteredIC(i1) = masterIC(i1);	% Resets the first looping IC to default
	end
end
	

%% Method 3: TKO for Nodes 4 & 8 -- Solver
% Solves when nodes 4, 8, and some other node gets knocked out
if nodeMethod == 3
	fixEqu(1) = 4;		% Next 4 lines will fix Node 4 & 8 to 0
	fixEqu(2) = 8;
	alteredIC(fixEqu(1)) = 0;
	alteredIC(fixEqu(2)) = 0;
	for i1 = 1:maxNode
		if i1== fixEqu(1)
			continue
		elseif i1 == fixEqu(2)
			continue
		else
			fixEqu(3) = i1;		% States equation to fix
			alteredIC(i1) = 0;	% Fixes the node's IC to 0
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			yData(intCount,:) = Y(end,:);	% Saves data
			xData(intCount,:) = [fixEqu(1), fixEqu(2), fixEqu(3)];
			alteredIC(i1) = masterIC(i1);	% Resets the IC to default values
			intCount = intCount + 1;
		end
	end
end


%% Method 4: Scaling Master Nodes -- Solver
% Solves when the Master Nodes are fixed
if nodeMethod == 4
	for i1 = 1:masterNode
		fixEqu(i1,:) = masterNodes(i1);		% This sets the master nodes to be fixed
		fprintf('\nRegarding Node %g\n',fixEqu(i1,:));
		multNode = input('Multiply said Node by = ');		% This sets the value of each master node
		alteredIC(fixEqu(i1,:)) = multNode*steadyValueP(fixEqu(i1,:));
	end
	
	[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
	
	yData(1,:) = Y(end,:);	% Saves data
	xData(1,:) = fixEqu;
	intCount = intCount + 1;
end


%% Method 5: Fix Nodes -- Solver
% This lets the user change the node to a desired value
if nodeMethod == 5;
	fix = input('How many nodes should be fixed? (0 -> 10) = ');
		if fix >= 1
			fixEqu = zeros(1,fix);
			for i1=1:1:fix
				fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
				fixEqu(i1) = fixNode(i1);		% Saves Node # to a matrix for use in the Function
				fixValue = input('Knockout with the Value = ');		% How the node should be modified
				alteredIC(fixNode(i1)) = fixValue;		% Replaces the value of Node # to the declared value
			end
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		else
		end	
end
    

%% Method 6: Sweep 1 Node -- Solver
% Solves a Custom Sweep of 1 Node
if nodeMethod == 6
	fixNode = input('Knockout Node Number = ');		% Node # to modify
		fixEqu(1) = fixNode;		% Saves Node # to a matrix for use in the Function
	maxValue = input('Select the max value your node should be multiplied by = ');
		avgValue = steadyValueP(fixNode);	% Retrieves the WT Value
		nodeLimit = maxValue * avgValue;	% Multiplies WT Value by Scaling Constant
		
	for i1 = 0:0.1:nodeLimit
		alteredIC(fixNode) = i1;		% Fixes value of Node #
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		yData(intCount,:) = Y(end,:);
		xData(intCount,1) = i1;
		intCount = intCount + 1;		
	end
end


%% Method 7: Sweep 2 Nodes -- Solver
% Solves a Custom Sweep of 2 Nodes
if nodeMethod == 7
	maxValue = input('Select the max value the two nodes should be multiplied by = ');
	for i1=1:1:2
		fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
		fixEqu(i1) = fixNode(i1);		% Saves Node # to a matrix for use in the Function
		nodeLimit(i1) = steadyValueP(fixNode(i1)) * maxValue;	% Multiplies WT Value by Scaling Constant
	end
		
	for i1 = 0:1:nodeLimit(1)
		alteredIC(fixNode(1)) = i1;		% Fixes value of Node 1
		fprintf('Node %g is at %g \n',fixNode(1),i1);
		for i2 = 0:1:nodeLimit(2)
			alteredIC(fixNode(2)) = i2;		% Fixes value of Node 2
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			Data(i1+1,i2+1,:) = Y(end,:);
			intCount = intCount + 1;		
		end
	end
end


%% Method 8: Manipulate Slave Nodes for Master Nodes -- Setup
% Creates the Matrices needed to solve for Method 8

if nodeMethod == 8
% Creates a Data Matrix of when each Master Node is knocked out
	for i1 = 1:masterNode
		alteredIC(masterNodes(i1)) = 0;		% Fixes the node's IC to 0
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		dataMatrix(i1,:) = Y(end,:);	% Saves SKO Values of Master Nodes in Data Matrix
		alteredIC = masterIC;
	end
	dataMatrix(masterNode+1,:) = steadyValueP(1,:);		% Saves WT Values in the last row of dataMatrix

% Creates a WT Matrix with WT Values
	for i1 = 1:masterNode+1
		wtMatrix(i1,:) = steadyValueP(:,1)';	% Creates a wtMatrix as large as dataMatrix
	end
	
% Creates a Difference Matrix by the difference of DataMatrix - WTMatrix
	diffMatrix = dataMatrix - wtMatrix;
	diffMatrix(end,:) = [];		% Eliminates the last row of zeros
			
% Creates Impact Matrix to reflect Identity Matrices
	impactMatrix = blkdiag(diffMatrix,-diffMatrix');
	impactMatrix((masterNode+1):end,1:maxNode)...
		= eye(maxNode);   %identity matrix for lower lft blk
	

%% Method 8: Manipulate Slave Nodes for Master Nodes -- Solver
% Solves by creating many Weight Matrices or using a User-Defined Matrix

% Applies Weight Matrix
	weightMatrix = weightValues(ones(masterNode,1),:);
	impactMatrix(1:masterNode,1:maxNode) = weightMatrix.*...
		impactMatrix(1:masterNode,1:maxNode);
	
% Prepares Vector for Left Hand Side of Equation
	dotProduct = diffMatrix.*desiredValues(ones(masterNode,1),:).*weightMatrix;
	vectorValues = sum(dotProduct,masterNode);
	vectorValues = [vectorValues;steadyValueP];
	
% Solves the System
	yData = linsolve(impactMatrix,vectorValues);	% The closest set of nodes to the actual 
	xData = yData(1:maxNode,:) - desiredValues';	% The diff between actual & desired values
end


%% Method 1: SKO for All -- Data Presentation
% Presents Data for Method 1
if nodeMethod == 1;
fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
	for i = 1:intCount
	fprintf('KO %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,1),yData(i,:));
	end
fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 2: DKO for All -- Data Presentation
% Presents Data for Method 2
if nodeMethod == 2;
fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
	for i = 1:intCount
	fprintf('KO %g & %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,:),yData(i,:));
	end
fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 3: TKO for Nodes 4 & 8 -- Data Presentation
% Presents Data for Method 3
if nodeMethod == 3;
fprintf('\nResults from yAxis Matrix for Node 4 & 8 \n');
	for i = 1:intCount
	fprintf('KO %g & %g & %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,:),yData(i,:));
	end
fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 4: Scaling Master Nodes -- Data Presentation
% Presents Data for Method 4
if nodeMethod == 4;
fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
	fprintf('Results: \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',yData(1,:));
	fprintf('Wildtype: \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 5: Fix Nodes -- Data Presentation
% Provides options for presenting data for Method 5
if nodeMethod == 5;
	
% This plots Time vs. Organism & Organism vs. Organism
	graph = input('\nShould plots be generated? (y/n) = ','s');
	if graph == 'y';
		plot(T,Y(:,1));
	else
	end
	
% This lets the user retrieve the Steady State Values
	steady = input('\nShould steady-state values be printed? (y/n) = ','s');
	if steady=='y'
		fprintf('The steadyState Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%f, ',Y(end,i));
			fprintf('/n/n');
		end
	else
	end
    
% This lets the user retrieve the Steady State Values
	steady = input('\nShould constants & alpha values be printed? (y/n) = ','s');
	if steady=='y'
		fprintf('The Constant Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%f, ',constantR(i));
			fprintf('/n');
		end
		fprintf('\n\nThe Alpha Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%f, ',newAlpha(1,i));
			fprintf('/n/n');
		end
	end
end


%% Method 6: Sweep 1 Node -- Data Presentation
% Presents Data for Method 6
if nodeMethod == 6;
nodeName = ['Node'];
	for i1 = 1:1:maxNode
		hold all;
		if i1 == fixNode
			continue
		else
			nodeStr = num2str(i1);
			nodeStr = strcat(nodeName,nodeStr);
			plot(xData(:,1),yData(:,i1),'--','DisplayName',nodeStr);
			legend('-DynamicLegend');
		end
	end
end


%% Method 7: Sweep 2 Nodes -- Data Presentation
% Presents Data for Method 7
if nodeMethod == 7;
	grid on
	fprintf('\nGenerating Plots. Note Nodes %g & %g are +1 larger in these plots\n',fixNode);
	
	subDiv = ceil(maxNode^0.5);		% Determines width & height of Subplots
	
	for i1 = 1:1:maxNode
		subplot(subDiv,subDiv,i1)
		mesh(Data(:,:,fixNode(1)), Data(:,:,fixNode(2)), Data(:,:,i1))	% This plots the data
		xlabel(['Node ',num2str(fixNode(1))]); ylabel(['Node ',num2str(fixNode(2))]);
		title(['Node ',num2str(i1)]);
	end
end


%% Method 8: Manipulate Slave Nodes for Master Nodes -- Data Presentation
if nodeMethod == 8;
	fprintf('Realistic Nodes:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',yData(1:maxNode));
	fprintf('Theoratical Nodes:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',desiredValues(1:maxNode));
	fprintf('Which is a difference of:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',xData);
end
	

%% Finish
% Achievement Unlocked.
fprintf('\nThis model is complete\n\n');