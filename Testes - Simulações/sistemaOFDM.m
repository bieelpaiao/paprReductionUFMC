s = rng(211);       % configura��o do gerador de n�meros aleat�rios

numFFT = 512;        % número de pontos para a FFT
subbandSize = 20;    % tamanho da sub-banda
numSubbands = 10;    % numSubbands*subbandSize <= numFFT

subbandOffset = (numFFT/2)-((subbandSize*numSubbands)/2); % centraliza as bandas no gr�fico do espectro

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 15;              % SNR in dB

%ajustes do gr�fico OFDM
hFig = figure;


% OFDM TRADICIONAL---------------------------------------------------------
% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada (matriz m x n, m = bitsPerSubCarrier*subbandSize, n = numSubbands)

%Loop over each subband
for bandIdx = 1:numSubbands

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1); %criando um sinal de entrada aleatório entre 0s e 1s
    % QAM Symbol mapper
    symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true); %criação do símbolo qam(sinal de entrada, ordem, argumentos do tipo nome=valor)
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
    
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize;
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];

    ifftOut = ifft(ifftshift(symbolsInOFDM));  
end

% Compute peak-to-average-power ratio (PAPR)
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
paprOFDM = pm(ifftOut);
disp(['PAPR do OFDM Tradicional = ' num2str(paprOFDM) ' dB']);
plotCCDF(pm);
hold on;

% OFDM PR�-CODIFICADO -----------------------------------------------------
np_v = [5 10 15];
for i=1:3
    %configura��o dos par�metros para a matriz de pr�-codifica��o
    nofdm = subbandSize;
    np = np_v(i);
    L = nofdm + np;
    beta = np/nofdm;
    pu0 = zeros(L, 1);
    puv = zeros(L,nofdm);
    
    %cria��o de pu,0
    for u = 0:L-1
        if u <= np
            pu0(u+1) = ((-1)^u/sqrt(nofdm))*sin(pi*u/(2*np));
        elseif u > np && u <= nofdm 
            pu0(u+1) = ((-1)^u/sqrt(nofdm));
        elseif u > nofdm && u <= L-1 
            pu0(u+1) = ((-1)^u/sqrt(nofdm))*cos(pi*(u-nofdm)/(2*np));
        end
    
        for v = 0:nofdm-1
            puv(u+1, v+1) = pu0(u+1)*exp((-1j*2*pi*(v)*(u))/nofdm);
        end
    end
    
    inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada (matriz m x n, m = bitsPerSubCarrier*subbandSize, n = numSubbands)
    
    for bandIdx = 1:numSubbands
    
        bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1); %criando um sinal de entrada aleatório entre 0s e 1s
        % QAM Symbol mapper
        symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
        'UnitAveragePower', true); %criação do símbolo qam(sinal de entrada, ordem, argumentos do tipo nome=valor)
        inpData(:,bandIdx) = bitsIn; % log bits for comparison
        
        symbolsCoded = puv * symbolsIn;
        
        % Pack subband data into an OFDM symbol
        offset = subbandOffset+(bandIdx-1)*subbandSize;
        symbolsInOFDM = [zeros(offset,1); symbolsCoded; ...
                         zeros(numFFT-offset-subbandSize, 1)];
    
        ifftOut = ifft(ifftshift(symbolsInOFDM));   
    end
    
    % Compute peak-to-average-power ratio (PAPR)
    pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
    paprOFDM = pm(ifftOut);
    disp(['PAPR do OFDM Pr�-codificado (' num2str((np/nofdm)*100) '%) = ' num2str(paprOFDM) ' dB']);
    plotCCDF(pm);
end

legend('OFDM Tradicional', 'OFDM Pr�-codificado (25%)', 'OFDM Pr�-codificado (50%)', 'OFDM Pr�-codificado (75%)');
xlabel('Pot�ncia relativa (dB acima da pot�ncia m�dia)');
ylabel('Probabilidade (%)')
title('CCDF')
hold off;




