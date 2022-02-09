%% ECE 5390 - Practicum 3 - Simulation Using the Frequency and Time Domains
%  G.Davis
%  02/10/2022

clc; clear; close all;

%% Differential Equations for Systems
syms s
de = laplace((1/(s+2))*(3/((s+4)*(s+5))));
sys1 = tf([1], [1 2]);
sys2 = tf([3], conv([1 4], [1 5]));
tf = sys1 * sys2;


