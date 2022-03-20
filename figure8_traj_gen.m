%clear all; clc; close all;
function [xtraj_unit,initialconditon] = figure8_traj_gen(indexcolloison)
global dt N_sim_end


    time_end = N_sim_end*dt;
    %indexcolloison = randi([120,180],1);
    patrol_period = indexcolloison*dt;
    ttraj = 0:dt:time_end;
    %% Figure 8
    centerloc_figure8 = [0;0;0];
    tlist = linspace(0,patrol_period,17);
    tlist = tlist(1:length(tlist)-1);
    xlist = [0 0.5 1 1.5 2 1.5 1 0.5 0 -0.5 -1 -1.5 -2 -1.5 -1 -0.5];
    ylist = [0 1 1.5 1 0 -1 -1.5 -1 0 -1 -1.5 -1 0 1 1.5 1];
    zlist = zeros(1,length(tlist));
    dup = ceil(time_end/8);
    count = 1;
    while count < dup
        tlist = [tlist, tlist+patrol_period*count];
        xlist = [xlist, xlist];
        ylist = [ylist, ylist];
        zlist = [zlist, zlist];
        count = count + 1;
    end
     ncut = min(find(tlist>=time_end));
     tlist = tlist(1:ncut);
     xlist = xlist(1:ncut);
     ylist = ylist(1:ncut);
     zlist = zlist(1:ncut);
    xlistq = interp1(tlist,xlist,ttraj);
    ylistq = interp1(tlist,ylist,ttraj);
    zlistq = interp1(tlist,zlist,ttraj);
    xtraj_unit = [xlistq;ylistq;zlistq]';
    initialconditon = [randn,0,randn,0,randn,0,0,0]';
end
% figure();
% plot3(xlist,ylist,zlist);
% figure();
% plot3(xlistq,ylistq,zlistq);
% 
% 
% figure();
% subplot(3,1,1);
% plot(tlist,xlist);
% hold on;
% plot(ttraj,xlistq);
% subplot(3,1,2);
% plot(tlist,ylist);
% hold on;
% plot(ttraj,ylistq);
% subplot(3,1,3);
% plot(tlist,zlist);
% hold on;
% plot(ttraj,zlistq);

