%clear all; clc; close all;
function [xtraj_unit,initialconditon] = hover_waypoint_gen(indexcolloison)
global dt N_sim_end
    time_end = N_sim_end*dt;
    ttraj = 0:dt:time_end;
    ztraj_unit = linspace(-3,1,length(ttraj));
    xtraj_unit = [[0;0]*ttraj;ztraj_unit];
    initialconditon = [xtraj_unit(1);0;xtraj_unit(2);0;xtraj_unit(3);0;0;0];
end
