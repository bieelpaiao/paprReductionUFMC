function pim = matrixGeneration(N_OFDM, Np_OFDM)
    %%
    N = N_OFDM;
    Np = Np_OFDM;
    L = N + Np;
    
    %%
    pi0 = complex(zeros(L, 1));
    pim = complex(zeros(L, N));
    
    %%
    i = 0:Np-1;
    pi0(1:Np) = (((-1).^i)/sqrt(N)).*sin((pi*i)/(2*Np));
    
    i = Np:N-1;
    pi0(Np+1:N) = ((-1).^i)/sqrt(N);
    
    i = N:L-1;
    pi0(N+1:end) = (((-1).^i)/sqrt(N)).*cos((pi*(i-N))/(2*Np));
    
    %%
    for i = 0:L-1
        for m = 0:N-1
            pim(i+1, m+1) = pi0(i+1)*exp((-1j*2*pi*(i)*(m))/N);
        end
    end
end