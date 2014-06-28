%% AMatrix Function: Introduction
% This MATLAB function calculates the A Matrix.
%		> The matrix elements describe the interactions between the nodes.
%		> Input Requirements: 1) inputMatrix & 2) N
%
% 1) inputMatrix: Matrix with normal circuit & perturbed data for 
%       all nodes.
%       Here, PER#1 refers to node values when node #1 is fixed to a value.
%             NOR refers to values of normal, unperturbed, network
%		This matrix should be in the form of:
%					node1	node2	node3	...	N
%		PER#1	
%		PER#2
%		PER#3
%		...
%		PER#N
%		NOR												 
%
% 2) N: The # of Nodes in this system.


%% AMatrix Function: Intitilization
function [A] = AMatrix(inputMatrix,N)
	A = zeros(N);	% Pre-Allocates an N x N Matrix for 'A'	

	
%% Creates a N+1 by N Matrix of WT Values
	P = inputMatrix(end,:);  % Obtains NOR from the last row of inputMatrix
	P = P(ones(N+1,1),:);	% Creates Matrix


%% Calculates A Matrix: Concept
% 1) The following section of code begins by duplicating
%		the NOR & Node Values for each iteration of the forLoop.
%		> This forLoop runs for each node that gets knocked out.
%
% 2) During kth iteration, we move P(k) to the left of the equation so 
%    P(k) = X-wtVal. 
%		>  This is done because the matrix equation which creates
%			a vector of A components {ex. A(1,1:N)} is 0 = A*(X - norVal).
%		>  Because this is a Homogenous System we'd get A = 0.
%		>  We can avoid this solution by conducting this step #2 and...	
%
% 3) We set norVal column 'k' to zero. This ensures A(k,k) = 1.
%
% 4) We remove the kth equation (row) because it's invalid.
%		>  This is because the 'k' node was fixed. Thus, it doesn't evolve 
%           according to the equation.
%
% 5) We solve for a new Right Side Matrix of (X - norVal).
%		>  Note: The NOR(last) row of this matrix is zero, except for X(k).
%		>  The subtraction cancels out most of these values because,
%			NOR = inputMatrix. It doesn't cancel X(k) because its value
%			still exists from Step#2 when P(k) was moved to the left.
%		>  The NOR Equation is now: P(k) = 0*A(k,j/=k) + X(k)*A(k,k).
%           P(k) = A(k,k) X(k). X(k) = P(k), so A(k,k) = 1.
%
% 6) Solves for the kth row of A Matrix using linsolve.

	
%% Calculates A Matrix: Code
for(k = 1:N)
	
	% 1) Duplicates the NOR & Node Values for the forLoop
	X = inputMatrix;
	norVal = P;

	% 2) Moves P(k) to the left hand side of the equation
	leftVal = P(:,k);  
	
	% 3) Sets norVal Column 'k' to zero
	norVal(:,k) = 0;

	% 4) Remove kth equation because it's not valid
	X(k,:) = [];
	norVal(k,:) = [];
	leftVal(k) = [];
	  
	% 5) Creates a new Right Side Matrix
	rightVal = X - norVal; % Subtracts NOR from Node Values except for X(k)
					  
	% 6) Solves for A Matrix 
	A(k,:) = linsolve(rightVal,leftVal);
	  
end