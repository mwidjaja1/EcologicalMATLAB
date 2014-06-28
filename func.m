%% Matthew Widjaja.
% Research Function Set.

% Instructions: Run these equations from program.m

function dy = func(t,y)
global fixEqu 

%% Initial Variables 
% This determines the characteristics of each organism

growthRate = [1; 2; 2; 1; 1; 2; 2; 1; 1; 2]; 
carryCap = [15; 15; 15; 15; 15; 15; 15; 15; 15; 15];  
compEffect = [0, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.30, 1.45, 0.30;... 
    0.75, 0, 0.90, 1.05, 1.30, 1.45, 0.30, 0.45, 0.60, 0.75;...
    1.30, 1.45, 0, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.30;...
    0.45, 0.60, 0.75, 0, 0.90, 1.05, 1.30, 1.45, 0.30, 0.45;...
    0.90, 1.05, 1.30, 1.45, 0, 0.30, 0.45, 0.60, 0.75, 0.90;...
    1.45, 0.30, 0.45, 0.60, 0.75, 0, 0.90, 1.05, 1.30, 1.45;...
    0.60, 0.75, 0.90, 1.05, 1.30, 1.45, 0, 0.30, 0.45, 0.60;...
    1.05, 1.30, 1.45, 0.30, 0.45, 0.60, 0.75, 0, 0.90, 1.05;...
    0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.30, 1.45, 0, 0.30;...
    0.75, 0.90, 1.05, 1.30, 1.45, 0.30, 0.45, 0.60, 0.75, 0];  
compEffect = compEffect./1.5;
a = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];   % Summation Alpha-Matrix

%% Calculation of Alpha Values
% This calculates the Alpha Value of each organism

for i=1:10  % Creates Alpha-Values for each organism
    for c=1:10
        if i == c;  % Skips adding Alpha if it affects itself
            a(i) = a(i);
        else % Adds Alpha for other organisms
            a(i) = a(i) + (compEffect(i,c) * y(c));
        end
    end
end

%% Generation of Equations
% This creates the most appropriate Diff-Equ for each organism

dy = zeros(10,1);
for i=1:10
    if any(i == fixEqu);    % If working on the Fixed Equation, we set it to a no-growth rate
        dy(i) = 0;
    else % Otherwise, we create the standard Equation
        dy(i) = (growthRate(i)* y(i)) * ((1/carryCap(i)) * (carryCap(i) - y(i) - (a(i))));
    end
end
end