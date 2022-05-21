function [ref_traj,xtraj,initialconditon,indexcolloison] = fun_const_spd_MC()
%%
global dt N_sim_end
ttraj = 0:dt:N_sim_end*dt;

%random initial conditions
xobsode0 = [rand*10,rand*10,rand*10]';
vobsode0 = [rand*5,rand*5,rand*5]';
mobs=rand+0.5;


xobstraj = xobsode0;
for i = 1:N_sim_end
    tspan = [ttraj(i), ttraj(i)+dt];
    [tobsode,xobsode]=ode45(@(t,y)fun_const_spd(t,vobsode0),tspan,xobsode0);
    xobsode0 = xobsode(end,:)';
    xobstraj = [xobstraj,xobsode0];
end
xtraj = xobstraj(1:3,:)';

indexcolloison = randi([140,180],1);
[reftraj_unit,initialconditon] =  hover_waypoint_gen(indexcolloison);%figure8_traj_gen(indexcolloison);
figure8centerloc = xtraj(indexcolloison,:);
ref_traj = figure8centerloc+reftraj_unit';%*(rand*2+0.5);
ref_traj = [ref_traj,zeros(N_sim_end+1,1)];
initialconditon = initialconditon + [figure8centerloc(1);0;figure8centerloc(2);0;figure8centerloc(3);0;0;0];
%
% figure();
% plot3(xtraj(:,1),xtraj(:,2),xtraj(:,3));
% hold on;
% plot3(ref_traj(:,1),ref_traj(:,2),ref_traj(:,3));
% hold on;
% plot3(initialconditon(1),initialconditon(3),initialconditon(5),'o');
end
