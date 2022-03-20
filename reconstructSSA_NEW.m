function G = reconstructSSA_NEW(Z,L,T)
    Zshifted = [];
        Zcut = Z;
        Zshift = NaN(L,T);
        for i = 1:L
            Zshift(i,i:T-L+i) = Zcut(i,:);
        end
        Zshifted = [Zshifted, Zshift];
    
    for i = 1:length(Zshifted(1,:))
        recon_states(1,i) = mean(rmmissing(Zshifted(:,i)));
    end
    G = recon_states(1,1:T);
 
end