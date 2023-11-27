%% Initialization
clear all, close all, clc

pos_init = [-3 -2 -1];
pos_final = [5 3 2];

tf = 8;
time_step = 0.01;

horizon = tf / time_step;

%% Quadrotor DDP

% Think we'll generally have pretty high Q_tf and Q Matrices, low R Matrix.

%% Quadcopter Recursive Model Predictive Control (rMPC)

%% Obstacle Avoidance with Discrete Barrier States