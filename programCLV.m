%% Matthew Widjaja.
% Research Program Runner.
% Instructions: This file & funcCLV.m must be in the same folder in order to be used.

clear all
format short

global maxNode, global alteredIC, global steadyValP,
global oldAlpha, global fixEqu, global newAlpha, global constantR
global pred_a, global effcient_e


%% -- Step 0: Clear Workspace Prompt --------------------------------------
% Prompts the user if the workspace & command window should be cleared.
% -------------------------------------------------------------------------
userQuery = input('Should Command Window be cleared? (y/n) = ','s');
if userQuery == 'y'
	clc
end


%% -- Step 1: Declare General Variables & Parameters ----------------------
% This function solves the ODEs with an RK45 method. 
%
% These globally defined variables are defined in this section:
%   maxNode     User Set    The amount of nodes/organisms present
%   alteredIC   Automated   This will hold the modified initial condtions
%                           of each node. When a node is fixed, the value
%                           said node is fixed to will be stored in this
%                           matrix.
%   steadyValP  User Set    Defines the steady state values of each node, 
%                           that is, the populations for which each node
%                           should ideally converge towards.
%   oldAlpha    User Set    These are the old alpha values which consists
%                           of several parameters and before predation is
%                           accounted for.
%
% These globally defined variables are defined in function file for the
%   Competitive Lotka-Volterra Model: fixEqu, newAlpha, constantR
%
% These globally defined variables are defined in function file for the
%   Competitive Lotka-Volterra Model: pred_a
%
% These variables are defined in this section are are used in other steps:
%   equName     User Set    The .m file (without the .m) containing the
%                           Lotka-Volterra Function which will be used.
%   intCount    Automated   This counts how many times iterations the model
%                           has ran, in order to build the results matrix.
%   maxTime     User Set    This defines how long the model will run for.
%
% These variables are defined in this section:
%   masterNode  User Set    This defines how quantity of master nodes.
%   masterNodes User Set    This matrix defines which nodes are master.
%   masterIC    Automated   This sets an initial condition of '1' for each
%                           node. This is used to reset alteredIC (which is
%                           modified when a node is fixed) after the model
%                           is completed in fixing a specific node.
% -------------------------------------------------------------------------

% Step 1A: State number of Nodes & Master Nodes
maxNode = 8;		% Amount of Nodes present
masterNode = 2;		% Amount of Master Nodes present
	masterNodes = [4 8];	% States which nodes are Master -- For Method 8
    
% Step 1B: State the duration of the model
maxTime = 7500;		% Max amount of time to use
equName = @funcNLV;	% Name of file w. Equations

% Step 1C: Initial Conditions
masterIC = ones(1,maxNode);     % Set master initial conditions
alteredIC = masterIC;           % Creates a matrix for the altered ICs

% Step 1D: Interger Counter
intCount = 0;

% Step 1E: Defines the variables for the competitive LV
steadyValP = [37.6583262579766, 46.6720836475590, 26.8674882660994, ... 
    26.0160062116315, 46.7387831740059, 44.3536550187196, ...
    26.5644981471383, 49.7072030311665];	
oldAlpha = [-0.0139, -0.0063, -0.0088, -0.0224, 0, 0, 0, -0.00637366326755053;
    -0.0025, -0.0189, -0.0169, -0.0122, 0, 0, 0, -0.00258912672726849; 
    -0.0011, -0.0133, -0.0274, -0.0201, 0, 0, 0, -0.00340351541750547; 
    -0.0076, -0.0056 ,-0.0144, -0.0315, 0, 0, 0,-0.0118562222381121; 
    0, 0, 0, -0.000789018955068564, -0.0364, -0.0067, -0.0033, -0.0078; 
    0, 0, 0, -0.00449527435703973 ,-0.007, -0.0349, -0.0036, -0.0341;  
    0, 0, 0, -0.00100095417983515, -0.0039, -0.0036, -0.0254, -0.0038; 
    0, 0, 0, -0.00361539381179086, -0.0044, -0.0091, -0.0084, -0.0141];	

% Step 1F: Defines the variables for the Niche LV model
pred_a = importdata('pred_a.txt');
effcient_e = 0.1;

% Step 1G: These parameters are used for Method 8
weightValue = 0.5;
for i1 = 1:maxNode
	weightValues(i1) = weightValue;
end
% weightValues = [0, 0, 1, 0, 1, 0, 0, 0];	
desiredValues = [18.9308, 17.9217, 0, 57.0000, 45.7847, 12.9408, 0.0010, 10.0000]; 


%% -- Step 2: Wildtype Data -----------------------------------------------
% This runs the model once, before any additional node modifications,
% to obtain the expected wildtype populations of each organism.
% -------------------------------------------------------------------------

fixEqu = 0;             % Fixes no equation
[T,W] = ode45(equName, [0 maxTime], masterIC );     % Solves the model
WT(1,:) = W(end,:);		% Saves WT data to matrix: WT
    

%% -- Step 3: Node Modification Prompt ------------------------------------
% The user is promoted to select how to modify each node (if any), 
% while running the model.
% -------------------------------------------------------------------------
fprintf('\nPossible Methods to run the Model:');
fprintf('\n0. Manual');
fprintf('\n1. Automated SKO');
fprintf('\n2. Automated DKO');
fprintf('\n3. Automated TKO');
fprintf('\n4. Scaled Master Nodes'); 
fprintf('\n5. Single Node Sweep');
fprintf('\n6. Double Node Sweep');
fprintf('\n7. Manipulate Slave Nodes for a Master Node');

fprintf('\n\n');
nodeMethod = input('Select a Method (0-7) = ');



%% -- QUERY D.4.0: Manual -------------------------------------------------
% This query lets the user run the ODE Function & fix specified nodes.
% 
% These variables will be identically used as previously documented: 
%   fixEqu, maxNode, constantR, newAlpha, maxTime
%
% The following parameters & variables will be modified in Method 0:
%   fixNode		User Set	States which nodes should be fixed.
%   alteredIC   User Set   Fixes the given node with the specified value. 
% -------------------------------------------------------------------------

if nodeMethod == 0;

% 0A: Lets the user change any number of nodes to a desired value
fix = input('How many nodes should be fixed? (0 -> 10) = ');
    if fix >= 1
        fixEqu = zeros(1,fix);
        for i1=1:1:fix
            fixNode(i1) = input('Knockout Node Number = ');	
            fixEqu(i1) = fixNode(i1);
            alteredIC(fixNode(i1)) = input('Knockout with the Value = ');
        end
    end
    
% 0B: Solves the model (using 'equName')
[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
	
% 0C: Lets the user display a line plot for each population
graph = input('\Display plots? (y/n) = ','s');
if graph == 'y';
    nodeName = ['Node'];
    
    % This automates the creation & naming of the line for each node
    for i1 = 1:1:maxNode
        hold all;
        nodeStr = num2str(i1);  
        nodeStr = strcat(nodeName,nodeStr);
        plot(T(:,1),Y(:,i1),'--','DisplayName',nodeStr);
        legend('-DynamicLegend');
    end
else
end
	
% 0D: Prints the Steady State, Constants, and new Alpha values.
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
    
    
%% -- Method 1: Automated SKO ---------------------------------------------
% This query does an automated Single Knock Out Routine.
%
% These globally defined variables will be used: fixEqu, equName, maxTime
%
% The following parameters & variables will be modified in Method 1:
%   alteredIC       Automatic	Fixes every node assigned to it as '0'
%   results         Results     Saves the population data
%
% The Algorithm used for Method 2 is:
%   1A. Fixes a node to 0 in an iterative loop from Node 0 to 'maxNode'
%       (maxNode is the last node)
%   1B. Solves the model (using the function file saved as 'equName').
%   1C. Saves the last timestep of data into results and the initial
%       conditions for the fixed node is reset to the original settings.
%   1D. The results are printed out to the terminal.
% -------------------------------------------------------------------------

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


%% -- Method 2: Automated DKO ---------------------------------------------
% This query does an automated Double Knock Out Routine.
%
% These variables will be identically used as previously documented: 
%   fixEqu, intCount, equName, maxTime
%
% The following parameters & variables will be modified in Method 2:
%   alteredIC       Automatic	Fixes every node assigned to it as '0'
%   results         Results     Saves the population data
%
% The Algorithm used for Method 2 is:
%   2A. Fixes a node to 0 in an iterative loop from Node 0 to 'maxNode-1'
%       (maxNode is the last node, and is included in the 2nd loop below).
%   2B. Fixes a second node to 0 in an iterative loop starting from the
%       node immediately following the node selected in the first loop in 
%       2A until the maxNode (ie. the last node).
%   2C. Solves the model (using the function file saved as 'equName').
%   2D. Saves the last timestep of data into results & the initial
%       conditions for 2B's fixed node is reset to the original settings.
%   2E. The results are printed out to the terminal.
%   2F. Resets the initial conditions of 2A's fixed node to the original
%       settings.
% -------------------------------------------------------------------------

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



%% -- Method 3: Automated TKO ---------------------------------------------
% This query does an automated Triple Knock Out Routine.
%
% These variables will be identically used as previously documented: 
%   fixEqu, intCount, equName, maxTime
%
% The following parameters & variables will be modified in Method 3:
%   alteredIC       Automatic	Fixes every node assigned to it as '0'
%   results         Results     Saves the population data
%
% The Algorithm used for Method 3 is:
%   3A. Both master nodes are fixed to 0
%   3B. A third non-master node is fixed to 0, one at a time
%   3C. Solves the model (using the function file saved as 'equName')
%   3D. Saves the last timestep of data into results and the initial
%       conditions for node 3 is reset to the original settings
%   3E. The results are printed out to the terminal
% -------------------------------------------------------------------------

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


%% -- Method 4: Scaled Master Nodes ---------------------------------------
% Each master node is multipled by a certain value that the user describes.
% From there, we solve the model and print the results to the terminal.
%
% These variables will be identically used as previously documented:
%   fixEqu, equName, alteredIC, maxTime
%
% The following parameters & variables will be modified in Method 4:
%   multNode	User Set	States the value the node is multiplied with.
%   alteredIC   Automatic   Multiplies the wild-type value of each fixEqu
%                           master node with the multNode value.
%
% The Algorithm used for Query 4 is:
%   4A. Select each master node for modification
%   4B. Fixes each master node by multiplying its wildtype value by a
%       user-specified value (multNode) and saving it in alteredIC
%   4C. Solves the model (using the function file saved as 'equName')
%   4D. Prints out results to the terminal
% -------------------------------------------------------------------------
     
if nodeMethod == 4

for i1 = 1:masterNode
    % 4A: Select each master node via a loop
    fixEqu(i1,:) = masterNodes(i1);
    
    % 4B: Fixes each MN by multiplying it by a user-specified value.
    fprintf('\nRegarding Node %g\n',fixEqu(i1,:));
    multNode = input('Multiply said Node by = ');
    alteredIC(fixEqu(i1,:)) = multNode * steadyValP(fixEqu(i1,:));
end
	
% 4C: Solves the Model (using 'equName')	
[T,Y] = ode45(equName, [0 maxTime], alteredIC );

% 4D: Prints out results
fprintf('Fix:\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',Y(end,:));
fprintf('WT:\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\t %7f\n',WT(1,:));
end


%% -- Method 5: Automated One-Node Sweep ----------------------------------
% This query lets the user do an automated sweep of one node as its value
% changes from 0 to a max value in specified stepSize increments.
%
% These variables will be identically used as previously documented:
%   fixEqu, steadyValP, equName, maxNode, intCount, maxTime
%
% The following parameters & variables will be modified in Method 5:
%   maxValue    User Set    Prompts the user for how the multiplier value
%                           of which the fixed node's WT value should be
%                           multiplied by.
%   nodeLimit	Automatic	Multiplies the Wildtype value of the fixed node
%                           by maxValue
%   alteredIC   Automatic   Fixes the node to value within 0 & nodeLimit in
%                           increments of 0.1.
%   timeData    Results     Stores the values of each timestep
%   result      Results     Stores the values of each node's population
%
% The Algorithm used for Query 5 is:
%   5A. Prompts the user for which one node to fix & saves it in fixEqu.
%   5B. Prompts the user to give the maxValue, the scaling constant for how
%       the fixed node should be multiplied by. This value is multiplied
%       with the wildtype value of said node and is saved as nodeLimit.
%   5C. Iteratively fixes said node from 0 to nodeLimit in 0.1 increments.
%   5D. Solves the model (using the function file saved as 'equName')
%   5E. Saves the data to the timeData & results matrices.
%   5F. Produces a dynamic graph to plot every node.
% -------------------------------------------------------------------------

if nodeMethod == 5
	
% 5A: Prompts the user for which node to fix
fixEqu(1) = input('Knockout Node Number = ');		
        
% 5B: Prompts for maxValue (Scaling Constant) & obtain the nodeLimit value
maxValue = input('The max value your node should be multiplied by = ');
nodeLimit = maxValue * steadyValP(fixEqu(1));

% 5C: Iteratively fixes said node from 0 -> nodeLimit
for i1 = 0:0.1:nodeLimit
    intCount = intCount + 1;	
    alteredIC(fixEqu(1)) = i1;
    
    % 5D: Solves the Model (using 'equName')
    [T,Y] = ode45(equName, [0 maxTime], alteredIC );
    fprintf('Node at %f\n',i1);
    
    % 5E: Saves the current time & population to timeData & results
    timeData(intCount,:) = i1;
    results(intCount,1) = Y(end,:);
end % End iterative node fixing

% 5F: Produces a graph to plot every node
nodeName = ['Node'];
for i1 = 1:1:maxNode
    hold all;
    nodeStr = num2str(i1);
    nodeStr = strcat(nodeName,nodeStr);
    plot(results(:,1),timeData(:,i1),'--','DisplayName',nodeStr);
    legend('-DynamicLegend');
end % End of graph
end


%% -- Method 6: Automated Two-Node Sweep ----------------------------------
% This query lets the user do an automated sweep of two nodes as each of
% the nodes have its values approach from 0 to its specified max value
% (which could be different for both) in a specified stepSize increment.
%
% These variables will be identically used as previously documented: 
%   fixEqu, steadyValP, equName, maxNode, intCount, maxTime
%
% The following parameters & variables will be modified in Method 5:
%   nodeLimit	Automatic	Multiplies the Wildtype value of the fixed node
%                           by maxValue
%   alteredIC   Automatic   Fixes the node to value within 0 & nodeLimit in
%                           increments of 0.1.
%   timeData    Results     Stores the values of each timestep
%   result      Results     Stores the values of each node's population
%
%
% The Algorithm used for Query 6 is:
%   6A. Prompts the user to give the maxValue, the scaling constant for how
%       the fixed node should be multiplied by.
%   6B. Prompts the user on which two nodes to fix. Its Wildtype Value is
%       multiplied with the sclaing constant and is saved as nodeLimit.
%   6C. Iteratively fixes the first node from 0 to the nodeLimit of the 
%       first node in increments of 1.0.
%   6D. Iteratively fixes the second node from 0 to the nodeLimit of the 
%       second node in increments of 1.0.
%   6E. Solves the model (using the function file saved as 'equName')
%   6F. Saves the data to the results matrices.
%   6G. Produces a dynamic graph to plot every node.
% -------------------------------------------------------------------------

if nodeMethod == 6
	
% 6A: Prompts the user for both node's maxValue
maxValue = input('Select the max value the two nodes should be multiplied by = ');
	
% 6B: Prompts for which 2 nodes to fix & calculates its nodeLimit    
for i1=1:1:2
    fixEqu(i1) = input('Fix Node Number = ');		
    nodeLimit(i1) = steadyValP(fixEqu(i1)) * maxValue;
end

% 6C: Iteratively fixes the first node from 0 -> nodeLimit
for i1 = 0:1:nodeLimit(1)
    alteredIC(fixEqu(1)) = i1;		
    fprintf('Node %f is at %f \n',fixEqu(1),i1);
    
    % 6D: Iteratively fixes the second node from 0 -> nodeLimit
    for i2 = 0:1:nodeLimit(2)
        alteredIC(fixEqu(2)) = i2;
        
        % 6E: Solves the Model (using 'equName')
        [T,Y] = ode45(equName, [0 maxTime], alteredIC );
        
        % 6F: Saves the current time & population to timeData & results
        result(i1+1,i2+1,:) = Y(end,:);
        intCount = intCount + 1;		
    end
end

% 6G: Produces a graph to plot every node
	grid on
	fprintf('\nGenerating Plots. Note Nodes %g & %g are +1 larger in these plots\n',fixEqu);

	subDiv = ceil(maxNode^0.5);		% Determines width & height of Subplots

	for i1 = 1:1:maxNode
		subplot(subDiv,subDiv,i1)
		mesh(result(:,:,fixEqu(1)), result(:,:,fixEqu(2)), result(:,:,i1));
		xlabel(['Node ',num2str(fixEqu(1))]); ylabel(['Node ',num2str(fixEqu(2))]);
		title(['Node ',num2str(i1)]);
	end
end



%% Method 8: Manipulate Slave Nodes for Master Nodes
if nodeMethod == 7

% Creates a Data Matrix of when each Master Node is knocked out
	for i1 = 1:masterNode
		alteredIC(masterNodes(i1)) = 0;		% Fixes the node's IC to 0
		[T,Y] = ode45(equName, [0 maxTime], alteredIC );	% Solves the model
		dataMatrix(i1,:) = Y(end,:);	% Saves SKO Values of Master Nodes in Data Matrix
		alteredIC = masterIC;
	end
	dataMatrix(masterNode+1,:) = steadyValP(1,:);		% Saves WT Values in the last row of dataMatrix

% Creates a WT Matrix with WT Values
	for i1 = 1:masterNode+1
		wtMatrix(i1,:) = steadyValP(:,1)';	% Creates a wtMatrix as large as dataMatrix
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
	vectorValues = [vectorValues;steadyValP];
	
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