%% Matthew Widjaja.
% Research Program Runner.
% Instructions: This file & func_v5.m must be in the same folder in order to be used.

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
% General Parameters to help model efficiencyx
maxNode = 8;		% Amount of Nodes present
masterNode = 2;		% Amount of Master Nodes present
	masterNodes = [4 8];	% States which nodes are Masters -- For Method 8
maxTime = 7500;		% Max amount of time to use
equName = @func_v5;		% Name of file w. Equations
masterIC = ones(1,maxNode);		% Set master initial conditions
alteredIC = masterIC;		% Creates a matrix for the altered ICs
intCount = 0;		% Integer Counter for some methods


%% General Data
% General Data to run the model
steadyValueP = [37.6583262579766, 46.6720836475590, 26.8674882660994,  26.0160062116315, 46.7387831740059, 44.3536550187196,  26.5644981471383, 49.7072030311665];					
oldAlpha = [-0.0139, -0.0063, -0.0088, -0.0224, 0, 0, 0, -0.00637366326755053; -0.0025, -0.0189, -0.0169, -0.0122, 0, 0, 0, -0.00258912672726849;  -0.0011, -0.0133, -0.0274, -0.0201, 0, 0, 0, -0.00340351541750547;  -0.0076, -0.0056 ,-0.0144, -0.0315, 0, 0, 0,-0.0118562222381121;  0, 0, 0, -0.000789018955068564, -0.0364, -0.0067, -0.0033, -0.0078; 0, 0, 0, -0.00449527435703973 ,-0.007, -0.0349, -0.0036, -0.0341;  0, 0, 0, -0.00100095417983515, -0.0039, -0.0036, -0.0254, -0.0038; 0, 0, 0, -0.00361539381179086, -0.0044, -0.0091, -0.0084, -0.0141];				

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
	fixEqu = 0;
	[T,W] = ode45(equName, [0 maxTime], masterIC );		% Solves the model
	WT(1,:) = W(end,:);		% Saves WT data


%% Method 1: SKO for All
if nodeMethod == 1
	
% Solves when all nodes get single knocked out
	for i1 = 1:maxNode
		intCount = intCount + 1;
		fixEqu = i1;	% States equation to fix
		alteredIC(i1) = 0;	% Fixes the node's IC to 0
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		yData(intCount,:) = Y(end,:);	% Saves data
		xData(intCount,:) = i1;
		alteredIC = masterIC;	% Resets the IC to default values
	end
	
% Presents Data for Method 1
	fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
		for i = 1:intCount
		fprintf('KO %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,1),yData(i,:));
		end
	fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 2: DKO for All
if nodeMethod == 2
	
% Solves when all nodes get double knocked out
	for i1 = 1:maxNode-1
		fixEqu(1) = i1;		% States equation to fix
		alteredIC(i1) = 0;	% Fixes the node's IC to 0
		for i2 = (i1+1):maxNode
			intCount = intCount + 1;
			fixEqu(2) = i2;		% States equation to fix
			alteredIC(i2) = 0;	% Fixes the node's IC to 0
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			yData(intCount,:) = Y(end,:);	% Saves data
			xData(intCount,:) = [i1, i2];
			alteredIC(i2) = masterIC(i2);	% Resets the second looping IC to default
		end
		alteredIC(i1) = masterIC(i1);	% Resets the first looping IC to default
	end
	
% Presents Data for Method 2
	fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
		for i = 1:intCount
		fprintf('KO %g & %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,:),yData(i,:));
		end
	fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 3: TKO for Master Nodes -- Solver
if nodeMethod == 3
	
% Fix the Master Nodes to '0'
	fixEqu(1) = masterNodes(1);	
	fixEqu(2) = masterNodes(2);
	alteredIC(masterNodes(1)) = 0;
	alteredIC(masterNodes(2)) = 0;
	
% Fixes a third node (via loop) to 0 automatically & solves the model
	for i1 = 1:maxNode
		if i1== fixEqu(1)	% This & the next 3 lines will ignore Equ 4 & 8. 
			continue
		elseif i1 == fixEqu(2)
			continue
		else
			intCount = intCount + 1;
			fixEqu(3) = i1;		% States equation to fix
			alteredIC(i1) = 0;	% Fixes the node's IC to 0
			[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			yData(intCount,:) = Y(end,:);	% Saves data
			xData(intCount,:) = [fixEqu(1), fixEqu(2), fixEqu(3)];
			alteredIC(i1) = masterIC(i1);	% Resets the IC to default values
		end
	end

% Presents Data for Method 3
	fprintf('\nResults from yAxis Matrix for Master Nodes %d & %d \n', fixEqu(1), fixEqu(2));
	for i = 1:intCount
		fprintf('KO %g & %g & %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',xData(i,:),yData(i,:));
	end
	fprintf('Wildtype \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 4: Scaling Master Nodes
if nodeMethod == 4

% Sets the values of the Master Nodes
	for i1 = 1:masterNode
		fixEqu(i1,:) = masterNodes(i1);		% This sets the master nodes to be fixed
		fprintf('\nRegarding Node %g\n',fixEqu(i1,:));
		multNode = input('Multiply said Node by = ');		% This sets the value of each master node
		alteredIC(fixEqu(i1,:)) = multNode*steadyValueP(fixEqu(i1,:));
	end
	
% Solves the model	
	[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model

% Presents Data for Method 4
	fprintf('\nResults from yAxis Matrix where row KO_0 = WT Value: \n');
	fprintf('Results: \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',Y(end,:));
	fprintf('Wildtype: \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n',WT(1,:));
end


%% Method 5: Fix Nodes
if nodeMethod == 5;

% This lets the user change the node to a desired value
	fix = input('How many nodes should be fixed? (0 -> 10) = ');
		if fix >= 1
			fixEqu = zeros(1,fix);
			for i1=1:1:fix
				fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
				fixEqu(i1) = fixNode(i1);		% Saves Node # to a matrix for use in the Function
				fixValue = input('Knockout with the Value = ');		% How the node should be modified
				alteredIC(fixNode(i1)) = fixValue;		% Replaces the value of Node # to the declared value
			end
		end
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
	
% This plots Time vs. Organism & Organism vs. Organism
	graph = input('\Display plots? (y/n) = ','s');
	if graph == 'y';
		nodeName = ['Node'];
		for i1 = 1:1:maxNode
			hold all;
			if i1 == fixNode
				continue
			else
				nodeStr = num2str(i1);		% This & the next 3 lines automates each node's plot
				nodeStr = strcat(nodeName,nodeStr);
				plot(T(:,1),Y(:,i1),'--','DisplayName',nodeStr);
				legend('-DynamicLegend');
			end
		end
	else
	end
	
% This prints the Steady State, Constants, and new Alpha values.
	fprintf('\nThe steadyState Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%f, ',Y(end,i));
		end
	fprintf('\n\nThe Constant Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
				fprintf('%f, ',constantR(i));
		end
	fprintf('\n\nThe Alpha Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
				fprintf('%f, ',newAlpha(1,i));
		end
	end


%% Method 6: Sweep 1 Node
if nodeMethod == 6
	
% Selects the node that will be fixed & to which value (via multiplication)
	fixNode = input('Knockout Node Number = ');		% Node # to modify
		fixEqu(1) = fixNode;		% Saves Node # to a matrix for use in the Function
	maxValue = input('Select the max value your node should be multiplied by = ');
		avgValue = steadyValueP(fixNode);	% Retrieves the WT Value
		nodeLimit = maxValue * avgValue;	% Multiplies WT Value by Scaling Constant

% Sweeps the specified node & solves the model
	for i1 = 0:0.1:nodeLimit
		intCount = intCount + 1;	
		alteredIC(fixNode) = i1;		% Fixes value of Node #
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		yData(intCount,:) = Y(end,:);
		xData(intCount,1) = i1;
	end

% Presents Data for Method 6
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


%% Method 7: Sweep 2 Nodes -- Solver
if nodeMethod == 7
	
% Selects the two nodes to fix and to which value (via multiplication)
	maxValue = input('Select the max value the two nodes should be multiplied by = ');
	for i1=1:1:2
		fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
		fixEqu(i1) = fixNode(i1);		% Saves Node # to a matrix for use in the Function
		nodeLimit(i1) = steadyValueP(fixNode(i1)) * maxValue;	% Multiplies WT Value by Scaling Constant
	end
	
% Solves the model
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

% Presents Data for Method 7
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



%% Method 8: Manipulate Slave Nodes for Master Nodes
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
	
% Presents Data for Method 8
	fprintf('Realistic Nodes:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',yData(1:maxNode));
	fprintf('Theoratical Nodes:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',desiredValues(1:maxNode));
	fprintf('Which is a difference of:\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n',xData);
end
	

%% Finish
% Achievement Unlocked.
fprintf('\nThis model is complete\n\n');