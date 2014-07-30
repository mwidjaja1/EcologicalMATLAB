%% Matthew Widjaja.
% Research Program Runner.
% Instructions: This file & funcCLV.m must be in the same folder in order to be used.

clear all
format short

global equName, global alteredIC, global fixEqu, global maxNode,
global newAlpha, global steadyValueP, global oldAlpha, global constantR


%% -- Step 0: Clear Workspace Prompt --
% Prompts the user if the workspace & command window should be cleared.
userQuery = input('Should Command Window be cleared? (y/n) = ','s');
if userQuery == 'y'
	clc
end


%% -- Step 1: Declare General Variables & Parameters --
% These define the structure and settings of the model

% Step 1A: State number of Nodes & Master Nodes
maxNode = 8;		% Amount of Nodes present
masterNode = 2;		% Amount of Master Nodes present
	masterNodes = [4 8];	% States which nodes are Masters -- For Method 8
    
% Step 1B: State the duration of the model
maxTime = 7500;		% Max amount of time to use
equName = @funcCLV;	% Name of file w. Equations

% Step 1C: Initial Conditions
masterIC = ones(1,maxNode);     % Set master initial conditions
alteredIC = masterIC;           % Creates a matrix for the altered ICs

% Step 1D: Interger Counter?
intCount = 0;		% Integer Counter for some methods

% Step 1E: Defines the ideal (i.e. steady-state) population values
steadyValueP = [37.6583262579766, 46.6720836475590, 26.8674882660994, ... 
    26.0160062116315, 46.7387831740059, 44.3536550187196, ...
    26.5644981471383, 49.7072030311665];	

% Step 1F: Defines the oldAlpha values before predation is accounted for
oldAlpha = [-0.0139, -0.0063, -0.0088, -0.0224, 0, 0, 0, -0.00637366326755053;
    -0.0025, -0.0189, -0.0169, -0.0122, 0, 0, 0, -0.00258912672726849; 
    -0.0011, -0.0133, -0.0274, -0.0201, 0, 0, 0, -0.00340351541750547; 
    -0.0076, -0.0056 ,-0.0144, -0.0315, 0, 0, 0,-0.0118562222381121; 
    0, 0, 0, -0.000789018955068564, -0.0364, -0.0067, -0.0033, -0.0078; 
    0, 0, 0, -0.00449527435703973 ,-0.007, -0.0349, -0.0036, -0.0341;  
    0, 0, 0, -0.00100095417983515, -0.0039, -0.0036, -0.0254, -0.0038; 
    0, 0, 0, -0.00361539381179086, -0.0044, -0.0091, -0.0084, -0.0141];				

% Step 1G: These parameters are used for Method 8
weightValue = 0.5;
for i1 = 1:maxNode
	weightValues(i1) = weightValue;
end
% weightValues = [0, 0, 1, 0, 1, 0, 0, 0];	
desiredValues = [18.9308, 17.9217, 0, 57.0000, 45.7847, 12.9408, 0.0010, 10.0000]; 


%% -- Step 2: Wildtype Data --
% This runs the model once, before any additional node modifications,
% to obtain the expected wildtype populations of each organism.

	fixEqu = 0;             % Fixes no equation
	[T,W] = ode45(equName, [0 maxTime], masterIC );	% Solves the model
	WT(1,:) = W(end,:);		% Saves WT data to matrix: WT
    

%% -- Step 3: Node Modification Prompt --
% The user is promoted to select how to modify each node (if any), 
% while running the model.
fprintf('\nPossible Methods to run the Model:');
fprintf('\n0. Manual');
fprintf('\n1. Automated SKO');
fprintf('\n2. Automated DKO');
fprintf('\n3. Automated TKO');
fprintf('\n4. Scaled Master Nodes'); 
fprintf('\n6. Single Node Sweep');
fprintf('\n7. Double Node Sweep');
fprintf('\n8. Manipulate Slave Nodes for a Master Node');

fprintf('\n\n');
nodeMethod = input('Select a Method (0-8) = ');


%% -- Method 0: Manual Fix Nodes --
% This method lets the user fix any number of node to any value, as
% well as giving the option to present the data in a number of ways.
if nodeMethod == 0;

% 0A: Lets the user change any number of nodes to a desired value
    fix = input('How many nodes should be fixed? (0 -> 10) = ');
		if fix >= 1
			fixEqu = zeros(1,fix);
			for i1=1:1:fix
				fixNode(i1) = input('Knockout Node Number = ');		% Node # to modify
				fixEqu(i1) = fixNode(i1);           % Saves Node # to a matrix for use in the Function
				fixValue = input('Knockout with the Value = ');		% How the node should be modified
				alteredIC(fixNode(i1)) = fixValue;	% Replaces the value of Node # to the declared value
			end
		end
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
	
% 0B: Lets the user display a line plot for each population
	graph = input('\Display plots? (y/n) = ','s');
	if graph == 'y';
		nodeName = ['Node'];
		for i1 = 1:1:maxNode
			hold all;
            nodeStr = num2str(i1);  % This & the next 3 lines automates each node's plot
            nodeStr = strcat(nodeName,nodeStr);
            plot(T(:,1),Y(:,i1),'--','DisplayName',nodeStr);
            legend('-DynamicLegend');
		end
	else
	end
	
% This prints the Steady State, Constants, and new Alpha values.
	fprintf('\nThe steadyState Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%7f, ',Y(end,i));
		end
	fprintf('\n\nThe Constant Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%7f, ',constantR(i));
		end
	fprintf('\n\nThe Alpha Values in the order of Node 1 --> 8 are:\n');
		for i = 1:maxNode
			fprintf('%7f, ',newAlpha(1,i));
		end
	end
    
    
%% -- Method 1: Automated SKO --
% i1 increments from Node 1 to maxNode. We knock out one node at a time
% to '0' and solve the model. The resulting values are saved & then the
% initial conditions are reset to knock out the next node.
if nodeMethod == 1
    
for i1 = 1:maxNode
    % 1A: Fixes one node to '0'
    fixEqu = i1;            % States equation to fix
    alteredIC(i1) = 0;      % Fixes the node's IC to 0
        
    % 1B: Solves the Model (using 'equName')
    [T,Y] = ode45(equName, [0 maxTime], alteredIC );
        
    % 1C: Saves the last timestep of data & resets alteredIC
    results(i1,:) = Y(end,:);	% Saves data
    alteredIC = masterIC;       % Resets the IC to default values
	
    % 1D: Prints out results per each SKO Node
    fprintf('KO %g\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',fixEqu,results(i1,:));
end %End Node 1's KO

fprintf('WT\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',WT(1,:));
end


%% -- Method 2: Automated DKO --
% i1 increments from Node 1 to maxNode-1. We knock out one node at a time.
% Meanwhile, i2 increments from Node i1+1 (the current value of i1+1) to
% maxNode. With two nodes knocked out, we then solve the model. The data
% is then saved & we reset the initial conditions to knock the next node.
if nodeMethod == 2
	
for i1 = 1:maxNode-1
    % 2A: Fixes one node to '0'
    fixEqu(1) = i1;         % States equation to fix
    alteredIC(i1) = 0.0000;	% Fixes the node's IC to 0
	
    for i2 = (i1+1):maxNode
        %2B: Fixes a second node to '0'
        intCount = intCount + 1;
        fixEqu(2) = i2;		% States equation to fix
        alteredIC(i2) = 0;	% Fixes the node's IC to 0
        
        %2C: Solves the Model (using 'equName')
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
			
        % 2D: Saves the last timestep of data & resets Node 2's IC
        results(intCount,:) = Y(end,:);
        alteredIC(i2) = masterIC(i2);
		
        % 2E: Prints out results for each DKO operation    
        fprintf('KO %g, %g\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',fixEqu,results(intCount,:));
    end %End Node 2's KO
        
    % 2F: Resets Node 1's IC
	alteredIC(i1) = masterIC(i1);
end %End of Node 1's KO

fprintf('All WT\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',WT(1,:));
end


%% -- Method 3: Automated TKO --
% Both master nodes (as stated in the masterNodes matrix of Step 1) are
% set to 0 for the entirety of this method. Meanwhile, i1 increments from
% Node 1 to maxNode & fixes said node to '0', so long as it isn't a master
% node itself. Once that's done, we solve the model & save its data. We
% then reset the initial conditions before knocking the next node.
if nodeMethod == 3
	
% 3A: Fixes both master nodes to '0'
fixEqu(1) = masterNodes(1);	
fixEqu(2) = masterNodes(2);
alteredIC(masterNodes(1)) = 0;
alteredIC(masterNodes(2)) = 0;
	
for i1 = 1:maxNode
    if i1~=fixEqu(1) || i1~=fixEqu(2)
        % 3B: Fixes a third (non-master) node to 0
        intCount = intCount + 1;
        fixEqu(3) = i1;		% States equation to fix
        alteredIC(i1) = 0;	% Fixes the node's IC to 0
        
        % 3C: Solves the Model (using 'equName')
        [T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
        
        % 3D: Saves the last timestep of data & resets Node 3's IC
        results(intCount,:) = Y(end,:);
        alteredIC(i1) = masterIC(i1);

        % 3E: Prints out results for each TKO operation
        fprintf('KO %g, %g, %g\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',fixEqu,results(i1,:));
    end
end %End of the Slave Node's KO

fprintf('All Wildtype\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',WT(1,:));
end


%% -- Method 4: Scaled Master Nodes --
if nodeMethod == 4

for i1 = 1:masterNode
    % 4A: Modifies each master node via a loop
    fixEqu(i1,:) = masterNodes(i1);
    
    % 4B: Fixes each MN by multiplying it by a user-specified value.
    fprintf('\nRegarding Node %g\n',fixEqu(i1,:));
    multNode = input('Multiply said Node by = ');
    alteredIC(fixEqu(i1,:)) = multNode*steadyValueP(fixEqu(i1,:));
end
	
% 4C: Solves the Model (using 'equName')	
[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model

% 4D: Prints out results
fprintf('Results:\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n', Y(end,:));
fprintf('Wildtype:\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',WT(1,:));
end


%% Method 6: Sweep 1 Node
if nodeMethod == 6
	
% 6A: Prompts the user for which node to fix
fixNode = input('Knockout Node Number = ');		
fixEqu(1) = fixNode;    % Saves Node # as the node to fix
        
% 6B: Fixes said node by multiplying it by a user-specified value
maxValue = input('Select the max value your node should be multiplied by = ');
avgValue = steadyValueP(fixNode);	% Retrieves the WT Value
nodeLimit = maxValue * avgValue;	% Multiplies WT Value by Scaling Constant

for i1 = 0:0.1:nodeLimit
    % 6C: Iteratively fixes said node from 0 -> nodeLimit in an 0.1 stepsize
    intCount = intCount + 1;	
    alteredIC(fixNode) = i1;
    
    % 6D: Solves the Model (using 'equName')
    [T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
    
    % 6E: Saves the current time & population to timeData & results
    timeData(intCount,:) = Y(end,:);
    results(intCount,1) = i1;
end % End iterative node fixing

% 6F: Produces a graph to plot every node
nodeName = ['Node'];
for i1 = 1:1:maxNode
    hold all;
    nodeStr = num2str(i1);
    nodeStr = strcat(nodeName,nodeStr);
    plot(results(:,1),timeData(:,i1),'--','DisplayName',nodeStr);
    legend('-DynamicLegend');
end % End of graph
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