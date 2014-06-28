%% Matthew Widjaja.
% Research Function Set.
% Instructions: Run program_v2.m to utilize these equations

function dy = func(t,y)
global fixEqu 
global maxNode
global constant
global newAlpha

%% Initial Variables 
% This determines the characteristics of each organism

steadyValueP = [0.0060, 0.0026, 0.0065, 0.0069, 0.0075, 0.0045, 0.0008, 0.0023]; 
oldAlpha = [-3.0289, -0.3701, -0.717, -0.448, 0, 0, 0, -0.21633505016943;...
    -2.0538, -1.43, -0.8654, -0.7504, 0, 0, 0, -0.193866116564711;...
    -0.0321, -0.3477, -1.0458, -1.3652, 0, 0, 0, -0.0398132096260681;...
    -0.7317, -0.1804, -0.9557, -1.7452, 0, 0 ,0 ,-0.228368899861282;...
    0, 0, 0, -0.150931605147117, -1.9238,-0.3892,-0.6804,-0.3006;...
    0, 0, 0, -0.15553106814765,-0.0058, -1.4085, -0.8097, -0.4383;...
    0, 0, 0, -0.134422580879091, -0.0613, -0.5699, -0.8766, -0.0381;...
    0, 0, 0, -0.041028739217153, -0.8943, -0.2053,-0.6153,-1.0412];

%% Constants (p)
% This calculates the Constants using alpha & steadyValue

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