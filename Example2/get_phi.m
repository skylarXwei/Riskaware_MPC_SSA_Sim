function Phi = get_phi(U,r,L)
    Phi = zeros(L-1,1);
    for i = 1:r
            Phi = Phi + U(1:(L-1),i)*U(L,i)';
    end
    v = norm(U(L,1:r));
    Phi = Phi/(1-v^2);   
    Phi = flipdim(Phi,1);
end