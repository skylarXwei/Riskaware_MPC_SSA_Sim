function [Urlist,rlist,Philist] = bootstrap_model(xtraj_training,H,L,Ntraining,rmax)
    covH = H*H';
    [U,~] = eig(covH);
    U = flip(U,2);
    V = [];
    dValue = [];
    V = H.'*U(:,1);
    Z = U(:,1)*V(:,1)';
    G = reconstructSSA_NEW(Z,L,Ntraining);
    value = norm(G - xtraj_training,2);
    i = 2;
    rcount = 1;
    while i<L && rcount<=rmax
        V = [V,H.'*U(:,i)];
        Z = U(:,1:i)*V(:,1:i)';
        G = reconstructSSA_NEW(Z,L,Ntraining);
        value(i) = norm(G - xtraj_training,2);
        dvalue(i-1) = abs((value(i) - value(i-1))/value(i));
        if dvalue(i-1) < 0.2
            rlist(rcount) = i;
            Urlist{rcount} = U(:,1:i);
            Philist(rcount,:) = get_phi(U,i,L);
            rcount = rcount+1;
        end
            i = i+1;
    end   
end