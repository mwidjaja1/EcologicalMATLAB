clear all; clc;

%% --- Matthew Widjaja: Niche Food Webs ---
% Purpose: This MATLAB code creates a Niche food web
%
% Credit:
% Matthew Widjaja:  Code Developer
% Jason Shulman:    Researcher for this model
% Russ Manson:      Reseracher for the LV Model
% Michael Widjaja:  Reminds me why I should have never grown up
%
%
%% -- Algorithm --
%   1. Obtain the number of species_S & create a niche_n value for each.
%   2. Create a beta-distributed random number (rand_x) for each species,
%       based on a constant beta value which is derived from the qty of
%       species and links.
%   3. Calculate the range_r values (with respect to niche_n & rand_x) and
%       the range's midpoint value (mid_c) for each species.
%   4. Compares each species(j) with species(i). If species(j) posesses a
%       niche_n value which is in the range of species(i)'s range (ie. in
%       the range of midRange < niche_n(j) < maxRange), then we say that
%       species(j) consumes species(i) & say pred_a(i,j)=1. Otherwise,
%       species(j) doesn't consume species(i), & pred_a(i,j) = 0.
%
%
%% --- List of Parameters & Variables ---
% Programming Parameters:
%   1. species_S = Qty of Species (from Step 1)
%       A user-prompted value for the qty of species
%   2. links_L = Qty of Links (from Step 2A)
%       A user-prompted value for the qty of links between species
%
% Model Variables:
%   1. niche_n[S] = Some random Niche Value (from Step 1)
%       A (uniform distribution) random value from 0-1.
%	2. range_r[S] = Niche Range (from Step 3)
%       The range of organisms who could be predated by organism i.
%	3. mid_c[S] = Range Midpoint (from Step 3)
%       Midpoint value of Ri. This is randomly generated from the range of
%       niche_n(i)/2 to niche(n).
%   4. pred_a[S,S] = Predation Values (from Step 4)
%       If organism 'j' will consume organism 'i', we'll get a[i,j] = 1.
%       If not, we'll get (a[i,j] = 0)
% 
% Beta Distribtion Parameters/Variables:
%   1. dirConn_C = Direct Connectance Value (from Step 2A)
%       Represents how connected an organism is to others with respect to
%       the number of links & species in the food web.
%   2. beta_B = Beta Value (from Step 2B)
%       This value was derived from 2*dirConn_c = (1/1+beta_B)
%   3. rand_x[S] = Beta-Distribution Values (from Step 2C)
%
% Additional Range Variables:
%   1. halfRange_r[S] = range_r[S]/2 (from Step 3)
%       Half of the range_r value which is used in further calculations.
%   2. minRange = Minimum range for species(i) (from Step 4)
%       Calculated from mid_c(i) - halfRange_r(i)
%   3. maxRange = Maximum range for species(i) (from Step 4)
%       Calculated from mid_c(i) + halfRange_r(i)


%% -- Step 1: Choose Ni --
% This step prompts the user for how many species (species_S) should be
% made, and then assigns a random niche value (niche_n) from 0-1 for each
% species.

species_S = input('Number of Species to Generate = ');
niche_n = rand(1,species_S);


%% -- Step 2: Beta Distribution Generator --
% This step declares the values needed to build an optimal Beta
% Distribution to generate random numbers for Step 3.

% 2A: dirConn_C = links_L/(species_S)^2
links_L = input('Number of Links to Generate = ');
dirConn_C = links_L/(species_S)^2;

% 2B: Solve for B if 2C=(1/1+B)
beta_B = (1/(2*dirConn_C))-1;

% 2C: Generate (beta-distributed) Random Numbers for each species
rand_x = betarnd(1,beta_B,1,species_S);


%% -- Step 3: Calculate Ri & Ci --
% This step calcuates the niche range (range_r) for each species and then
% determines each range's midpoint (mid_c). To do this, we will need to
% multiply the Niche Value of each species by each species' random beta-
% distributed number which was derived in Step 2. We also note that the
% node with the smallest niche_n value will always have a range_r of 0.

% 3A: Allocates Arrays
range_r = zeros(1,species_S);
halfRange_r = zeros(1,species_S);
mid_c = zeros(1,species_S);

for i=1:1:species_S
    % 3B: If node i has low niche value, it's a basal species w. r(i) = 0
    if min(niche_n) == niche_n(i)
        range_r(i) = 0;
        basal_b = i;
    
    % 3C: Otherwise, r(i) = n & some beta-distributed random value
    else
        range_r(i) = niche_n(i) * rand_x(i);
    end
    
    % 3D: Calculating the mid_c value of each node
    halfRange_r(i) = range_r(i)/2;
    mid_c(i) = (niche_n(i)-halfRange_r(i))*rand(1) + halfRange_r(i);
end


%% -- Step 4: Calculate Aij --
% This step calcuates if organism 'j' will consume organism 'i' (a(i,j)=1)
% or not (a(i,j)=0), by determining if organism j's niche_n value is in the
% range as defined by minRange & maxRange.

pred_a = zeros(species_S,species_S);

for j=1:1:species_S
    minRange = mid_c(j) - halfRange_r(j);
    maxRange = mid_c(j) + halfRange_r(j);
    for i=1:species_S
        %if i==j
        %    pred_a(i,j) = 0;
        if minRange <= niche_n(i) <= maxRange
            fprintf('Yes: Range of %g-%g with niche %g\n', minRange,maxRange,niche_n(i));
            pred_a(i,j) = 1;
            %pred_a(i,j) = 0;
        else
            fprintf('No: Range of %g-%g with niche %g\n', minRange,maxRange,niche_n(i));
            pred_a(i,j) = 0;
            %pred_a(i,j) = 1;
        end
    end 
end

%% -- Step 5: Save Aij & Basal Species --
% pred_a is saved to a text file.
save('niche','pred_a','basal_b')
