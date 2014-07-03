clear all; clc;

%% --- Matthew Widjaja: Cascade Food Webs ---
% Purpose: This MATLAB code creates a cascade food web
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
%       species(i) consumes species(j) & say pred_a(i,j)=1. Otherwise,
%       species(i) doesn't consume species(j), & pred_a(i,j) = 0.
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
%       A (uniform distribution) random value from 0-1 which indicates
%       organism i's likelihood of being consumed (where 1 = least likely)
%	2. range_r[S] = Niche Range (from Step 3)
%       The range of organisms who could be predated by organism i.
%	3. mid_c[S] = Range Midpoint (from Step 3)
%       Midpoint value of Ri. This is randomly generated from the range of
%       niche_n(i)/2 to niche(n).
%   4. pred_a[S,S] = Predation Values (from Step 4)
%       The rate on if organism 'i' will consume organism 'j' (so that
%       a[i,j] = 1) or not (a[i,j] = 0)
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

% Step 2A: dirConn_C = links_L/(species_S)^2
links_L = input('Number of Links to Generate = ');
dirConn_C = links_L/(species_S)^2;

% Step 2B: Solve for B if 2C=(1/1+B)
beta_B = (1/(2*dirConn_C))-1;

% Step 2C: Generate (beta-distributed) Random Numbers for each species
rand_x = betarnd(1,beta_B,1,species_S);


%% -- Step 3: Calculate Ri & Ci --
% This step calcuates the niche range (range_r) for each species and then
% determines each range's midpoint (mid_c). To do this, we will need to
% multiply the Niche Value of each species by each species' random beta-
% distributed number which was derived in Step 2.

range_r = ones(1,species_S);
halfRange_r = ones(1,species_S);
mid_c = ones(1,species_S);

for i=1:1:species_S
    range_r(i) = niche_n(i) * rand_x(i);
    halfRange_r(i) = range_r(i)/2;
    mid_c(i) = (niche_n(i)-halfRange_r(i))*rand(1) + halfRange_r(i);
end


%% -- Step 4: Calculate Aij --
% This step calcuates if organism 'i' will consume organism 'j' (a(i,j)=1)
% or not (a(i,j)=0), by determining if organism j's niche_n value is in the
% range (range_r+mid_

pred_a = ones(species_S,species_S);

for i=1:1:species_S
    minRange = mid_c(i) - halfRange_r(i);
    maxRange = mid_c(i) + halfRange_r(i);
    for j=1:1:species_S
        if i==j
            pred_a(i,j) = 0;
        elseif minRange <= niche_n(j) <= maxRange
            fprintf('Yes: Range of %g-%g with niche %g\n', minRange,maxRange,niche_n(j));
            pred_a(i,j) = 1;
        else
            fprintf('No: Range of %g-%g with niche %g\n', minRange,maxRange,niche_n(j));
            pred_a(i,j) = 0;
        end
    end 
end
