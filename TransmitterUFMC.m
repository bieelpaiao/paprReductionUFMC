% Definição dos Parâmetros Iniciais ---------------------------------------

s = rng(211);   % gerador de números aleatórios (manter repetitividade)

numFFT = 512;   % número de pontos da FFT
subbandSize = 20;   % tamanho de cada sub-banda
numSubbands = 10;   % número de sub-bandas
subbandOffset = (numFFT/2)-((subbandSize*numSubbands)/2); % offset para centralizar o espectro
filterLen = L + nofdm - 1;    % tamanho do filtro Chebyshev
slobeAtten = 40;    % atenuação do lóbulo lateral, em dB
bitsPerSubCarrier = 4;   % bits por subportadora

% Inicialização dos vetores e do filtro -----------------------------------
prototypeFilter = chebwin(filterLen, slobeAtten); %filtro com janela de Dolph-Chebyshev
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada
txSig = complex(zeros(numFFT+filterLen-1, 1)); %sinal da transmissão

% Definição dos parâmetros para a matriz de pré-codificação ---------------
nofdm = subbandSize;    % número de subportadoras por símbolo OFDM
np = subbandSize/2;    % subportadoras extras (sobrecarga)
L = nofdm + np;     % linhas da matriz de pré-codificação
pu0 = zeros(L, 1);    % índice p(u, 0)
puv = zeros(L,nofdm);    % matriz completa

% Criação da matriz de pré-codificação
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

% Processamento dos dados para cada sub-banda -----------------------------
for bandIdx = 1:numSubbands

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);    % bits de entrada
    symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', 'UnitAveragePower', true); % modulação dos bits de entrada
    
    symbolsCoded = puv * symbolsIn;    % pré-codificação
    
    % Juntando cada sub-banda para formar o símbolo OFDM
    offset = subbandOffset+(bandIdx-1)*subbandSize;    % deslocamento no espectro
    symbolsInOFDM = [zeros(offset,1); symbolsIn; zeros(numFFT-offset-subbandSize, 1)];    % criação do símbolo OFDM

    ifftOut = ifft(ifftshift(symbolsInOFDM));    % aplicação da IFFT para transmissão do sinal no domínio do tempo -> x[n]

    % Filtragem de cada sub-banda
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT*((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );     % janela para filtragem -> h[n]  
    filterOut = conv(bandFilter,ifftOut);    % saída do filtro -> convolução entre x[n] e h[n]

    txSig = txSig + filterOut;    % sinal a ser transmitido     
end

% Registro do PAPR do sinal UFMC pré-codificado
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);    % objeto pm para medição do PAPR
paprUFMC = pm(txSig);    % registro do PAPR do sinal UFMC
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);    % mostrar na tela

% Registro do PAPR do sinal OFDM pré-codificado
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);     % objeto pm para medição do PAPR
paprOFDM = pm(ifftOut);    % registro do PAPR do sinal OFDM
disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);    % mostrar na tela