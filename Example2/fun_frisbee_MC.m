function [ref_traj,xtraj,initialconditon,indexcolloison]  = fun_frisbee_MC()
%%
global dt N_sim_end
%%  Frisbee flight dynamics as per Sara Hummel's PhD thesis 
format short 
global mfrisbee gfrisbee Ia Id A d rho  
global CLo CLa CDo CDa CMo CMa CRr           
global CMq CRp CNr 
%%  define non-aerodynamic parameters (based on Sara Hummel's paper)
mfrisbee = 0.075;   % Kg 
gfrisbee = 9.7935;  % m/s^2 
A = 0.1963;   % m^2 
d = 2*sqrt(A/pi);  % diameter 
rho = 1.23; % Kg/m^3 
Ia  = 0.002352;  % moment of inertia about the spinning axis 
Id  = 0.001219; % moment of inertia about the planar axis'

%%  THE THREE ESTIMATED COEFFICIENTS
%CMq= -0.005,     CRp =-0.0055,   CNr =  0.0000071       % short (three) flights 
 CMq= -1.44E-02 ; CRp =-1.25E-02; CNr = -3.41E-05;       % long flight f2302 

%%  THE seven COEFFICIENTS estimated from three flights   
CLo=  0.3331; 
CLa=  1.9124; 
CDo=  0.5769; 
CDa=  1.85; 
CMo= -0.0821;   
CMa=  0.4338; 
CRr=  0.00171; 
%%  angle and angular velocities in rad and rad/sec respectively
%%  phi   = angle about the x axis  phidot       = angular velocity
%%  theta = angle about the y axis    thetadot    = angular velocity
%%  gamma = angle about the z axis   gd(gammadot)   = angular velocity          
%%  Common release conditions:  
%%    thetao = 0.192;  speedo = 14; betao = 0.105; gd=50  
%%  Define Simulation set initial conditions, enter your choosen values here: 
thetao =  -randi([10,40])*pi/180;%-.192;    % initial pitch angle 
speedo = 3+rand*5;%13.7;  % magnitude, m/sec 
betao  = randi([20,80])*pi/180;%.105;     % flight path angle in radians between velocity vector and horizontal 
gd= 300+randi([0,300]);            % initial spin 

alphao = thetao - betao;          % inital alpha 
vxo = speedo * cos(betao);        % velocity component in the nx direction
vzo = (speedo * sin(betao));     % velocity component in the nz direction  
                                 %    (note: nz is positive down) 
%%  x0= vector of initial conditions 
%%  x0= [   x    y     z  vx  vy  vz   phi theta   phidot  thetadot gd   gamma]  
%%  First set:
x0= [ rand*10 rand*10 rand*50  vxo 0 vzo  0   thetao  0.001   0.001    gd  0];
%%  Second set: 
tfinal = N_sim_end*dt;  % length of flight 
dtt = dt;
tspan=[0:dtt:tfinal];
options=[];
%%  Calls the ODE and integrate the frisbee equations of motions in the  
%%    subroutine, discfltEOM.m  
[t,x]=ode45(@discfltEOM,tspan,x0,options); 

xtraj = x(:,1:3);

indexcolloison = randi([140,180],1);
[reftraj_unit,initialconditon] =  hover_waypoint_gen(indexcolloison);%figure8_traj_gen(indexcolloison);
figure8centerloc = xtraj(indexcolloison,:);
ref_traj = figure8centerloc+reftraj_unit';%*(rand*2+0.5);
ref_traj = [ref_traj,zeros(N_sim_end+1,1)];
initialconditon = initialconditon + [figure8centerloc(1);0;figure8centerloc(2);0;figure8centerloc(3);0;0;0];

% indexcolloison = randi([140,200],1);
% [reftraj_unit,initialconditon] = figure8_traj_gen(indexcolloison);
% figure8centerloc = xtraj(indexcolloison,:);
% ref_traj = figure8centerloc+reftraj_unit;%*(rand*2+0.5);
% ref_traj = [ref_traj,zeros(N_sim_end+1,1)];
% initialconditon = initialconditon + [figure8centerloc(1);0;figure8centerloc(2);0;figure8centerloc(3);0;0;0];
%%
end
