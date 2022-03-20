function dqout = quad_dyn(t,qin,Uin,Ac,Bc)
    global m g
    dqout = Ac*qin + Bc*Uin;
    dqout(6) = dqout(6) - g;
end