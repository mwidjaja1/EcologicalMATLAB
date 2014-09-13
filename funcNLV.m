%% -- Matthew Widjaja: Niche LV Func --------------------------------------
% This code works in conjunction with programNLV.m and/or sortingHat.m in
% order to derive the diff equs needed to solve a Niche Lotka Volterra
% function. Run program.m first to use this function.

function dy = funcNLV(t,y)
global maxNode, global fixEqu
global pred_a, global effic_e, global basal_b 
global growth_b, global interact_a, global rand_a

%% -- Step 1: Calculates interact_a ---------------------------------------
% This uses pred_a from nicheWeb.m in order to create interact_a, which is
% a matrix describing the interactions between organisms.
% 
% Variables or Arrays used throughotu this code:
%   interact_a[] Automatic  This contains the combined growth & death rates
%                           for each organism. This is derived from the
%                           binary values of pred_a[]
%
% Variables or Arrays used in this step alone:
%   basal_b      Automatic  Created in nicheWeb.m, this single value notes
%                           which species had the lowest niche value, and
%                           thus, becomes the basal species.
%   pred_a[]	 Automatic	Created in nicheWeb.m, this states that if
%                           pred_a(i,j) = 1, then j eats i. Otherwise, if
%                           pred_a(i,j) = 0, then j does not eat i.
%   rand_v       Automatic  A random value from 0 to 1 to assist in the
%                           calculation rates of a cannibal species.
%   effic_e      User Set   Parameter defining the required amount of prey
%                           in order to produce a predator. Default = 0.1.
%   
% Algorithm for Step 1:
%   1A. An interact_a matrix is allocated with the size: maxNode x maxNode.
%   1B. We create an iteration loops for 1 > maxNode to check every node's
%       value in interact_a.
%       1C. If we are analyzing the basal species (as noted by basal_b from
%           nicheWeb.m), we set its interact_a(i,i) = -0.01 as stated in
%           the provided research paper.
%       1D. If we are analyzing a species that is a cannibal (as noted if
%           pred_a(i,i) = 1), then we modify its interact_a(i,i) value.
%       1E. Regardless of if 1C or 1D were met, we then create a loop to
%           compare said node (chosen in 1B's loop) to every other node.
%           >>> If pred_a(i,j) = 1, then we set interact_a(i,j) to some
%               random value between -1 & 0 and this indicates that j hurts
%               i. We then set interact_a(j,i) = effic_e * interact_a(i,j)
%               and this indicates that i helps j.
% -------------------------------------------------------------------------

% 1A: An interact_a matrix is allocated
interact_a = zeros(1,maxNode);

% 1B: We loop for each row in interact_a
for (i = 1:1:maxNode)
    
    % 1C: Sets interact_a(i,i) = -0.01 if Node i is a basal species
    if (i == basal_b);
        interact_a(i,i) = -0.01;  
    end %End of If Statement
    
    % 1D: Sets interact_a(i,i) if Node i is canabalistic: pred_a(i,i) = 1
    if (pred_a(i,i) == 1);
        rand_v = -rand(1);
        interact_a(i,i) = rand_v - (effic_e*rand_v);
    end %End of If Statement
        
    % 1E: Sets interact_a(i,j) & (j,i) if Node j eats i: pred_a(i,j)=1
    for (j=1:1:maxNode)
        if (pred_a(i,j) == 1 && i~=j);
                interact_a(i,j) = -rand(1);
                interact_a(j,i) = -effic_e * interact_a(i,j);
        end %End of If Statement
    end %End of 'j' Loop
        
end %End of 'i' Loop
   

%% -- Step 2: Calculates growth_b -----------------------------------------
% This gives each node a growth_b value as some random value from 0 to 1.
% -------------------------------------------------------------------------
growth_b = -rand(1,maxNode);
growth_b(basal_b) = rand(1);


%% -- Step 3: Generation of Equations -------------------------------------
% This creates the most appropriate Diff-Equ for each organism.

% Variables or Arrays used in this step:
%   sigma     Automated     The product of interact_a(i,j) & y(j) or x(j) 
%                           as we iterate over 'j' nodes from 0 to maxNode
%                           to solve for dy(i).
%
% Algorithm for Step 3:
%   3A. A matrix called dy is created to store the equation results
%   3B. We iterate through each node with as 'i' and reset 'sigma'
%   3C. In that loop, if the node we're working on is fixed, we set its
%       dy value to 0.
%   3D. In that loop, if the node we're working on isn't fixed, we solve
%       for its equation.
%       3D1. We iterate over each node again as 'j'. We do this because we
%            need to solve for ?a_ij*X_j over 'j' nodes where X & y equates
%            to population.
%       3D2. We then solve for dX_i/dt = X_i*(b_i+?a_ij*X_j)
% -------------------------------------------------------------------------

% 3A: A maxNode x 1 matrix is created to store the equation results 
dy = zeros(maxNode,1);

% 3B: Each node is iterated through
for(i = 1:maxNode)					   
    sigma = 0;      % Resets or Initiates Sigma
    
    % 3C: If the node we're working on is fixed, we set its value to 0
    if any(i == fixEqu);
        dy(i) = 0;
    
    % 3D: If the node we're working on isn't fixed, we solve its equation
    else
        %3D1: We iterate for each node again as j, summing the product of
        %     its interact_a with node 'i' & the population of node 'j'.
        for(j = 1:maxNode)
            sigma = sigma + (interact_a(i,j) * y(j));
        end % End of Node 'j' Loop
        
        %3D2: We solve for dX_i/dt = X_i * (b_i + ?a_ij * X_j)
        dy(i) = y(i) * (growth_b(i) + sigma);
    
    end %End of If Statement
end %End of Node 'i' Loop
end