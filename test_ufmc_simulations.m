clear all;
clc;

numSimulations = 1000;  % N�mero de simula��es para calcular a CCDF
PAPR_values = zeros(numSimulations, 1); % Array para armazenar os valores do PAPR

for sim = 1:numSimulations
    numFFT = 512;        % número de pontos para a FFT
    subbandSize = 20;    % tamanho da sub-banda
    numSubbands = 10;    % numSubbands*subbandSize <= numFFT

    %configuração dos parâmetros para a matriz de pré-codificação
    nofdm = subbandSize;
    np = subbandSize/2;
    L = nofdm + np;
    beta = np/nofdm;
    pu0 = zeros(L, 1);
    puv = zeros(L,nofdm);
    y = zeros(L, 1);

    subbandOffset = (numFFT/2)-(((subbandSize*numSubbands))/2); % centraliza as bandas no gráfico do espectro

    % Dolph-Chebyshev window design parameters
    filterLen = 43;      % similar to cyclic prefix length
    slobeAtten = 40;     % side-lobe attenuation, dB

    bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
    snrdB = 15;              % SNR in dB

    % Design window with specified attenuation
    prototypeFilter = chebwin(filterLen, slobeAtten);

    % Transmit-end processing
    %  Initialize arrays
    inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada (matriz m x n, m = bitsPerSubCarrier*subbandSize, n = numSubbands)
    txSig = complex(zeros(numFFT+filterLen-1+np, 1)); %sinal da transmissão (tamanho = numFFT+filterLen-1)

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
            puv(u+1, v+1) = pu0(u+1)*exp(-1i*2*pi*(v+1)*(u+1)/nofdm);
        end
    end

    %Loop over each subband
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

        % Filter for each subband is shifted in frequency
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
        filterOut = conv(bandFilter,ifftOut);

        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        txSig = txSig + filterOut;     
    end

    % Calcule o PAPR para o sinal UFMC atual
    power_peak = max(abs(txSig).^2); 
    power_average = mean(abs(txSig).^2);
    PAPR = 10 * log10(power_peak / power_average);
    
    PAPR_values(sim) = PAPR; % Armazena o valor do PAPR
    
    sim = sim + 1;
end

% Calcule a CCDF a partir dos valores do PAPR
PAPR_sorted = sort(PAPR_values, 'descend');
CCDF = (1:numSimulations) / numSimulations;

% Plote o gr�fico PAPR x CCDF
figure;
semilogy(PAPR_sorted, CCDF);
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF');
title('CCDF of UFMC PAPR');
