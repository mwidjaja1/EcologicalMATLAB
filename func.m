%% Matthew Widjaja.
% Research Function Set.
% Instructions: Run program_v3.m to utilize these equations

function dy = func(t,y)
global fixEqu 
global maxNode
global newAlpha
global steadyValueP


%% Initial Alpha Values 
% Species Competition Values (Unitless/Dimensionless)
oldAlpha = [-0.1636, -0.426, -0.1146, -0.1041, 0, 0, 0, -0.1723; ...
	-0.0393, -0.7424, -0.1717,	-0.0843, 0, 0, 0, -0.1416; ...
	-0.1751, -0.3305, -0.5004, -0.407, 0, 0 ,0, -0.2361; ...
	-0.1623, -0.0278, -0.2633, -0.4735, 0, 0, 0, -0.2517; ...
	0, 0, 0, -0.1601, -0.3113, -0.409, -0.3461, -0.3886; ...
	0, 0, 0, -0.1681, -0.1243, -0.5388, -0.2226, -0.22; ...
	0, 0, 0, -0.2088, -0.0128, -0.1183, -0.3811, -0.1343; ...
	0, 0, 0, -0.1394,-0.0897, -0.099, -0.0856, -0.3288];

%% Constants
% This calculates the Constants using oldAlpha & steadyValueP

constantR = zeros(1,maxNode);		% Generates Constants Matrix
for(i1 = 1:maxNode)						% Populates Constants Matrix
	for(i2 = 1:maxNode)
		newAlpha(i1,i2) = -oldAlpha(i1,i2) / steadyValueP(i1);
		constantR(i1) = constantR(i1) + (newAlpha(i1,i2) * steadyValueP(i2));
	end
end

%% Generation of Equations
% This creates the most appropriate Diff-Equ for each organism

dy = zeros(maxNode,1);	% Generates Equation Matrix
	for(i1 = 1:maxNode)					   
		if any(i1 == fixEqu);	% This fixes an Equation's growth if requested in program_v2.m
			dy(i1) = 0;
		else	% Otherwise, an Equation is made
			dy(i1) = constantR(i1) * y(i1);	% ...By multiplying the Constant by Y
			for(i2 = 1:maxNode)	% ...And subtracting Alpha
				dy(i1) = dy(i1) - (newAlpha(i1,i2) * y(i1) * y(i2));
			end
		end
	end
end