%% -- Matthew Widjaja: Niche Network Sorter -------------------------------
% This code creates niche networks to find a network with the steady state
% properties desired for our Ecological Network.

global maxNode, global fixEqu, global equName
global pred_a, global effic_e, global basal_b
global growth_b, global interact_a, global rand_a


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
% These variables are defined in this section are are used in other steps:
%   maxTime     User Set    This defines how long the model will run for.
%   maxNode     User Set    The amount of nodes/organisms present
%   equName     User Set    The .m file (without the .m) containing the
%                           Lotka-Volterra Function which will be used.
%
% These globally defined variables are generated in nicheWeb.m:
%   pred_a      Automated   This data is derived from the niche.mat which
%                           was made in nicheWeb.m. This defines the 
%                           interactions between 2 organisms in binary.
%   effic_e     User Set    This efficiency parameter is a value greater
%                           than 0 and defines the amount of prey needed
%                           to produce a predator.
%   basal_b     Automated   This data is derived from the niche.mat which
%                           was made in nicheWeb.m. This defines which
%                           species are basal species.
%
% On the niche LV, these variables are defined in funcNLV.m:
%   interact_a, growth_b, sigma
%
% These variables are defined in this section:
%   masterIC    Automated   This sets a random initial condition from 0-1
%                           for each node. 
% -------------------------------------------------------------------------
% 1A: Defines the variables used throughout the model
maxTime = 500;		% Max amount of time to use

% 1B: Defines the variables for the niche LV model
equName = @funcNLV;   % Name of file w. Equations
load('niche');        % The Niche Values as obtained from nicheWeb.m 
effic_e = 0.1;

% 1C: Initial Conditions
maxNode = length(pred_a);
masterIC = rand(1,maxNode);     % Set master initial conditions



%% -- Step 2: Wildtype Data -----------------------------------------------
% This runs the model once, to obtain the expected wildtype populations
% -------------------------------------------------------------------------
[T,W] = ode45(equName, [0 maxTime], masterIC );     % Solves the model


%% -- Step 3: Analyze Difference ------------------------------------------
% We analyze the results to see how many nodes aren't extinct
%
% These variables are used in this step alone:
%   passNode     Automatic      Counts how many nodes didn't go extinct
% -------------------------------------------------------------------------
WT(1,:) = W(end,:);                     % Saves WT data to matrix: WT
passNode = 0;

for i=1:maxNode
    if (0.001 < WT(1,i));
        passNode = passNode + 1;
    else
        continue;
    end %End of If Statement
fprintf('Nodes Passed: %g\n',passNode); 
end %End of For Statement


%% -- Step 4: Plot Results ------------------------------------------------
% This creates a plot to show the populations of every organism
% -------------------------------------------------------------------------
growth_b = growth_b;
interact_a = interact_a;
rand_a = rand_a;

nodeName = ['Node'];
for i1 = 1:1:maxNode
    hold all;
    nodeStr = num2str(i1);  
    nodeStr = strcat(nodeName,nodeStr);
    plot(T(:,1),W(:,i1),'--','DisplayName',nodeStr);
    legend('-DynamicLegend');
end %End of For Statement
