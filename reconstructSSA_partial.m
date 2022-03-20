function G = reconstructSSA_partial(Z,L)
    G = zeros(1,L);
    Zflip = flip(Z,2);
    for i = 1:L-1
        G(i) = mean(diag(Zflip,-i));
    end
end