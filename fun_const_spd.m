%Function to forward propagate the balll trajectory

function dy=fun_der(t,y)
%y=[x,y,z,xdot,ydot,zdot
dy=[y(1);
    y(2);
    y(3)];
end