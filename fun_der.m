%Function to forward propagate the balll trajectory

function dy=fun_der(t,y,Cd,m,d)
%y=[x,y,z,xdot,ydot,zdot
g=-9.81;
v=norm(y(4:6));
dy=[y(4);
    y(5);
    y(6);
    0-0.5*1.225*Cd*pi*d^2/(4*m)*v*y(4);
    0-0.5*1.225*Cd*pi*d^2/(4*m)*v*y(5);
    g-0.5*1.225*Cd*pi*d^2/(4*m)*v*y(6)];
end