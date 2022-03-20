%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  File:    discfltEOM.m
%%  By:      Sarah Hummel  
%%  Date:    July 2003 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot=discfltEOM(~,x) 
% Equations of Motion for the frisbee 
% The inertial frame, xyz = forward, right and down
global m g Ia Id A d rho  
global CLo CLa CDo CDa CMo CMa CRr           
global CMq CRp CNr 
% x = [ x y z vx vy vz f th  fd  thd  gd gamma]
%       1 2 3 4  5  6  7  8   9   10  11  12 
%% give states normal names 
vx = x(4); 
vy = x(5); 
vz = x(6); 
f  = x(7); 
th = x(8); 
st = sin(th); 
ct = cos(th); 
sf = sin(f); 
cf = cos(f); 
fd = x(9); 
thd= x(10);
gd  = x(11); 
     
%% Define transformation matrix 
%% [c]=[T_c_N] * [N]
T_c_N=[ct st*sf -st*cf; 0 cf sf; st -ct*sf ct*cf]; 
%% [d]=[T_d_N] * [N] 
%T_d_N(1,:)=[cg*ct  sg*cf+sf*st*cg  sf*sg-st*cf*cg];
%T_d_N(2,:)=[ -sg*ct cf*cg-sf*sg*st sf*cg+sg*st*cf];
%T_d_N(3,:)=[ st -sf*ct cf*ct]
     
c3=T_c_N(3,:);      % c3 expressed in N frame
     
%% calculate aerodynamic forces and moments 
%% every vector is expressed in the N frame
vel = [vx vy vz];   %expressed in N 
vmag = norm(vel); 
   
vc3=dot(vel,c3);     % velocity (scalar) in the c3 direction 
vp= [vel-vc3*c3];     % subtract the c3 velocity component to get the velocity vector  
% projected onto the plane of the disc, expressed in N 
alpha = atan(vc3/norm(vp)); 
Adp = A*rho*vmag*vmag/2; 
uvel  = vel/vmag;            % unit vector in vel direction, expressed in N 
uvp   = vp/norm(vp);      % unit vector in the projected velocity direction, expressed in N 
ulat  = cross(c3,uvp); % unit vec perp to v and d3 that points to right, right? 
%% first calc moments in uvp (roll), ulat(pitch) directions, then express in n1,n2,n3 
omegaD_N_inC = [fd*ct thd  fd*st+gd];       % expressed in c1,c2,c3 
omegaD_N_inN = T_c_N'*omegaD_N_inC';      % expressed in n1,n2,n3 
        
omegavp   = dot(omegaD_N_inN,uvp);         
omegalat  = dot(omegaD_N_inN,ulat);        
omegaspin = dot(omegaD_N_inN,c3);            % omegaspin = p1=fd*st+gd 
 
    CL  = CLo + CLa*alpha; 
    alphaeq = -CLo/CLa;    % this is angle of attack at zero lift 
    CD  = CDo + CDa*(alpha-alphaeq)*(alpha-alphaeq); 
    CM=CMo + CMa*alpha; 
     
    %CRr= CRr*d*omegaspinv/2./vmagv';
    %CRr= CRr*sqrt(d/g)*omegaspinv;    % this line produces NaN, so leave it in Mvp equation 
    %Mvp = Adp*d* (CRr*d*omegaspin/2/vmag  + CRp*omegavp)*uvp;   % expressed in N 
    Mvp = Adp*d* (sqrt(d/g)*CRr*omegaspin  + CRp*omegavp)*uvp;   % expressed in N 

  
lift  = CL*Adp; 
drag  = CD*Adp; 
ulift = -cross(uvel,ulat);          % ulift always has - d3 component 
udrag = -uvel;
Faero = lift*ulift + drag*udrag;     % aero force in N 
FgN   = [ 0 0 -m*g]';                 % gravity force in N 
F = Faero' + FgN; 
Mlat  = Adp*d* (CM + CMq*omegalat)*ulat;    % Pitch moment expressed in N 
Mspin = [0 0 +CNr*(omegaspin)];              % Spin Down moment expressed in C 
M = T_c_N*Mvp' + T_c_N*Mlat' + Mspin';     % Total moment expressed in C 
         
% set moments equal to zero if wanted... 
% M=[0 0 0];       
   
% calculate the derivatives of the states 
xdot = vel'; 
xdot(4)  = (F(1)/m);     %accx 
xdot(5)  = (F(2)/m);     %accy
xdot(6)  = (F(3)/m);     %accz 
xdot(7)  = fd;
xdot(8)  = thd;
xdot(9)  = (M(1) + Id*thd*fd*st - Ia*thd*(fd*st+gd) + Id*thd*fd*st)/Id/ct; 
xdot(10) = (M(2) + Ia*fd*ct*(fd*st +gd) - Id*fd*fd*ct*st)/Id; 
fdd=xdot(9); 
xdot(11) = (M(3) - Ia*(fdd*st + thd*fd*ct))/Ia; 
xdot(12) = x(11); 
