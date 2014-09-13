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
%   maxNode     User Set    The amount of nodes/organisms present
%   maxTime     User Set    This defines how long the model will run for.
%   intCount    Automated   This counts how many times iterations the model
%                           has ran, in order to build the results matrix.
%   masterNode  User Set    This defines how quantity of master nodes.
%   equName     User Set    The .m file (without the .m) containing the
%                           Lotka-Volterra Function which will be used.
%i
% If using Niche LV: These globally defined variables are defined:
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
% On the niche LV, these variables are defined in funcCLV.m:
%   interact_a, growth_b, sigma
%
% These variables are defined in this section:
%   masterNodes User Set    This matrix defines which nodes are master.
%   masterIC    Automated   This sets an initial condition of '1' for each
%                           node. This is used to reset alteredIC (which is
%                           modified when a node is fixed) after the model
%                           is completed in fixing a specific node.
% -------------------------------------------------------------------------
% 1A: Defines the variables used throughout the model
maxTime = 5000;		% Max amount of time to use
passNode = 0;
modelRun = 0;

% 1B: Initial Conditions
masterIC = ones(1,maxNode);     % Set master initial conditions

% 1C: Defines the variables for the niche LV model
maxNode = 50;
equName = @funcNLV;   % Name of file w. Equations
load('niche');        % Obtained from nicheWeb.m 
effic_e = 0.1;

% 1A: An interact_a matrix is allocated
interact_a = 500 * ones(1,maxNode);

% 1B: We loop for each row in interact_a
for(i = 1:maxNode)
    for(j = 1:maxNode) 
    
    % 1F: If pred_a(i,j)=1 & i?j, then interact_a(i,j) & (j,i) is set
    (pred_a(i,j) == 1 && i~=j);
        interact_a(i,j) = -rand(1);
        interact_a(j,i) = -effic_e * interact_a(i,j);
        
    % 1C: If Node i is a basal species, then interact_a(i,i) = -0.01
    if (i == basal_b);
        interact_a(i,i) = 0;    
    
    % 1D: If pred_a(i,i) = 1, then interact_a(i,i)
    elseif (pred_a(i,i) == 1);
        rand_v = -rand(1);
        interact_a(i,i) = rand_v - (effic_e*rand_v);
    	
    
    
    else
        continue
    end
    end
end