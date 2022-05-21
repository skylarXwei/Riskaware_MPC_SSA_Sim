function [Urlist,rlist,Philist] = bootstrap_model_with_randomQR(xtraj_training,H,L,Ntraining,rmax)
    %covH = H*H';
    %[U,~] = eig(covH);
    
    %U = flip(U,2);

for k = 1:L
   temp = rand(length(H(1,:)),1);
   Omega(:,k) = temp/norm(temp);
end
Y = H*Omega;
[U,~] = qr(Y,0);
H_tilt = U*(U'*H);
G = reconstructSSA_NEW(H_tilt,L,Ntraining);
value = norm(G - xtraj_training,2);
    i = 2;
    rcount = 1;
    while i<L && rcount<=rmax
        %V = [V,H.'*U(:,i)];
        %Z = U(:,1:i)*V(:,1:i)';
        H_tilt = U(:,1:i)*(U(:,1:i)'*H);
        G = reconstructSSA_NEW(H_tilt,L,Ntraining);
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


