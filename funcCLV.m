%% Matthew Widjaja.
% Research Function Set.
% Instructions: Run programCLV.m to utilize these equations

function dy = func_v5(t,y)
global maxNode, global fixEqu
global constantR,  global newAlpha, global oldAlpha, global steadyValP 

%% Constants
% This calculates the Constants using oldAlpha & steadyValueP

constantR = zeros(1,maxNode);		% Generates Constants Matrix
for(i1 = 1:maxNode)						% Populates Constants Matrix
	for(i2 = 1:maxNode)
		newAlpha(i1,i2) = -oldAlpha(i1,i2) / steadyValP(i1);
		constantR(i1) = constantR(i1) + (newAlpha(i1,i2) * steadyValP(i2));
	end
end

%% Generation of Equations
% This creates the most appropriate Diff-Equ for each organism

dy = zeros(maxNode,1);	% Generates Equation Matrix
	for(i1 = 1:maxNode)					   
		if any(i1 == fixEqu);	% This fixes an Equation's growth if requested in the Program
			dy(i1) = 0;
		else	% Otherwise, an Equation is made
			dy(i1) = constantR(i1) * y(i1);	% ...By multiplying the Constant by Y
			for(i2 = 1:maxNode)	% ...And subtracting Alpha
				dy(i1) = dy(i1) - (newAlpha(i1,i2) * y(i1) * y(i2));
			end
		end
	end
end