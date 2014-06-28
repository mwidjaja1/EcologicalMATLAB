%% Matthew Widjaja.
% Research Steady State Grapher

clear all
clc
format short

%% Data when Node 4 is Knocked-Out
% Obtained from knocking out Node 4 at the specified values

LF1 = [0.004925, 0.002287, 0.015664, 0.000000, 0.007755, 0.004981, 0.001548, 0.001816]; % = 0.000
LF2 = [0.005080, 0.002331, 0.014335, 0.001000, 0.007718, 0.004912, 0.001440, 0.001886]; % = 0.001
LF3 = [0.005236, 0.002377, 0.013007, 0.002000, 0.007681, 0.004842, 0.001331, 0.001956]; % = 0.002
LF4 = [0.005392, 0.002423, 0.011679, 0.003000, 0.007644, 0.004772, 0.001223, 0.002026]; % = 0.003
LF5 = [0.005548, 0.002468, 0.010351, 0.004000, 0.007607, 0.004702, 0.001115, 0.002096]; % = 0.004
LF6 = [0.005704, 0.002514, 0.009023, 0.005000, 0.007570, 0.004633, 0.001006, 0.002167]; % = 0.005
LF7 = [0.005860, 0.002559, 0.007695, 0.006000, 0.007533, 0.004563, 0.000898, 0.002237]; % = 0.006
LF8 = [0.006001, 0.002601, 0.006507, 0.006897, 0.007500, 0.004500, 0.000800, 0.002300]; % = WT

%% Modifiying the Data Set per Organism
% Creation of a For-Loop to Augment the Matrix

for i = 1:8
    Node4(:,i) = [LF1(i), LF2(i), LF3(i), LF4(i), LF5(i), LF6(i), LF7(i), LF8(i)];
end

%% Plotting of Data Set
xPlot = 3; yPlot = 3;
subplot(xPlot,yPlot,1)
    plot(Node4(:,4),Node4(:,1))
    title('Org 4 vs. Org 1')
subplot(xPlot,yPlot,2)
    plot(Node4(:,4),Node4(:,2))
    title('Org 4 vs. Org 2')
subplot(xPlot,yPlot,3)
    plot(Node4(:,4),Node4(:,3))
    title('Org 4 vs. Org 3')
subplot(xPlot,yPlot,4)
    plot(Node4(:,4),Node4(:,5))
    title('Org 4 vs. Org 4')
subplot(xPlot,yPlot,5)
    plot(Node4(:,4),Node4(:,6))
    title('Org 4 vs. Org 5')
subplot(xPlot,yPlot,6)
    plot(Node4(:,4),Node4(:,7))
    title('Org 4 vs. Org 7')
subplot(xPlot,yPlot,7)
    plot(Node4(:,4),Node4(:,8))
    title('Org 4 vs. Org 8')

