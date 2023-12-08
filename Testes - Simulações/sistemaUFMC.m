s = rng(211);       % configuração do gerador de números aleatórios

numFFT = 512;        % nÃºmero de pontos para a FFT
subbandSize = 20;    % tamanho da sub-banda
numSubbands = 10;    % numSubbands*subbandSize <= numFFT

subbandOffset = (numFFT/2)-((subbandSize*numSubbands)/2); % centraliza as bandas no gráfico do espectro

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 15;              % SNR in dB

% Design window with specified attenuation
prototypeFilter = chebwin(filterLen, slobeAtten);

%ajustes do gráfico UFMC
hFig = figure;

% UFMC TRADICIONAL---------------------------------------------------------
% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada (matriz m x n, m = bitsPerSubCarrier*subbandSize, n = numSubbands)
txSig = complex(zeros(numFFT+filterLen-1, 1)); %sinal da transmissão (tamanho = numFFT+filterLen-1)

%Loop over each subband
for bandIdx = 1:numSubbands

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1); %criando um sinal de entrada aleatÃ³rio entre 0s e 1s
    % QAM Symbol mapper
    symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true); %criaÃ§Ã£o do sÃ­mbolo qam(sinal de entrada, ordem, argumentos do tipo nome=valor)
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
    
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize;
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];

    ifftOut = ifft(ifftshift(symbolsInOFDM));

    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                 ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
    filterOut = conv(bandFilter,ifftOut);

    % Plot power spectral density (PSD) per subband
    % [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
    %                       numFFT*2, 1, 'centered'); 
    % plot(f,10*log10(psd)); 

    % Sum the filtered subband responses to form the aggregate transmit
    % signal
    txSig = txSig + filterOut;     
end

% Compute peak-to-average-power ratio (PAPR)
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
paprUFMC = pm(txSig);
disp(['PAPR do UFMC Tradicional = ' num2str(paprUFMC) ' dB']);
plotCCDF(pm);
hold on;


% set(hFig, 'Position', figposition([20 50 25 30]));
% hold off;


% UFMC PRÉ-CODIFICADO -----------------------------------------------------
np_v = [5 10 15];
for i=1:3
    %configuração dos parâmetros para a matriz de pré-codificação
    nofdm = subbandSize;
    np = np_v(i);
    L = nofdm + np;
    beta = np/nofdm;
    pu0 = zeros(L, 1);
    puv = zeros(L,nofdm);
    
    %criação de pu,0
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
    txSig = complex(zeros(numFFT+filterLen-1+np, 1)); %sinal da transmissão (tamanho = numFFT+filterLen-1)
    
    for bandIdx = 1:numSubbands
    
        bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1); %criando um sinal de entrada aleatÃ³rio entre 0s e 1s
        % QAM Symbol mapper
        symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
        'UnitAveragePower', true); %criaÃ§Ã£o do sÃ­mbolo qam(sinal de entrada, ordem, argumentos do tipo nome=valor)
        inpData(:,bandIdx) = bitsIn; % log bits for comparison
        
        symbolsCoded = puv * symbolsIn;
        
        % Pack subband data into an OFDM symbol
        offset = subbandOffset+(bandIdx-1)*subbandSize;
        symbolsInOFDM = [zeros(offset,1); symbolsCoded; ...
                         zeros(numFFT-offset-subbandSize, 1)];
    
        ifftOut = ifft(ifftshift(symbolsInOFDM));
    
        % Filter for each subband is shifted in frequency
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
        filterOut = conv(bandFilter,ifftOut);
    
        % Plot power spectral density (PSD) per subband
        % [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
        %                       numFFT*2, 1, 'centered'); 
        % plot(f,10*log10(psd)); 
    
        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        txSig = txSig + filterOut;     
    end
    
    % Compute peak-to-average-power ratio (PAPR)
    pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
    paprUFMC = pm(txSig);
    disp(['PAPR do UFMC Pré-codificado (' num2str((np/nofdm)*100) '%) = ' num2str(paprUFMC) ' dB']);
    plotCCDF(pm);
end

legend('UFMC Tradicional', 'UFMC Pré-codificado (25%)', 'UFMC Pré-codificado (50%)', 'UFMC Pré-codificado (75%)');
xlabel('Potência relativa (dB acima da potência média)');
ylabel('Probabilidade (%)')
title('CCDF')
hold off;




