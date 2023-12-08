% Defini��o dos Par�metros Iniciais ---------------------------------------

s = rng(211);   % gerador de n�meros aleat�rios (manter repetitividade)

numFFT = 512;   % n�mero de pontos da FFT
subbandSize = 20;   % tamanho de cada sub-banda
numSubbands = 10;   % n�mero de sub-bandas
subbandOffset = (numFFT/2)-((subbandSize*numSubbands)/2); % offset para centralizar o espectro
filterLen = L + nofdm - 1;    % tamanho do filtro Chebyshev
slobeAtten = 40;    % atenua��o do l�bulo lateral, em dB
bitsPerSubCarrier = 4;   % bits por subportadora

% Inicializa��o dos vetores e do filtro -----------------------------------
prototypeFilter = chebwin(filterLen, slobeAtten); %filtro com janela de Dolph-Chebyshev
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands); %dados de entrada
txSig = complex(zeros(numFFT+filterLen-1, 1)); %sinal da transmiss�o

% Defini��o dos par�metros para a matriz de pr�-codifica��o ---------------
nofdm = subbandSize;    % n�mero de subportadoras por s�mbolo OFDM
np = subbandSize/2;    % subportadoras extras (sobrecarga)
L = nofdm + np;     % linhas da matriz de pr�-codifica��o
pu0 = zeros(L, 1);    % �ndice p(u, 0)
puv = zeros(L,nofdm);    % matriz completa

% Cria��o da matriz de pr�-codifica��o
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
    symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', 'UnitAveragePower', true); % modula��o dos bits de entrada
    
    symbolsCoded = puv * symbolsIn;    % pr�-codifica��o
    
    % Juntando cada sub-banda para formar o s�mbolo OFDM
    offset = subbandOffset+(bandIdx-1)*subbandSize;    % deslocamento no espectro
    symbolsInOFDM = [zeros(offset,1); symbolsIn; zeros(numFFT-offset-subbandSize, 1)];    % cria��o do s�mbolo OFDM

    ifftOut = ifft(ifftshift(symbolsInOFDM));    % aplica��o da IFFT para transmiss�o do sinal no dom�nio do tempo -> x[n]

    % Filtragem de cada sub-banda
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT*((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );     % janela para filtragem -> h[n]  
    filterOut = conv(bandFilter,ifftOut);    % sa�da do filtro -> convolu��o entre x[n] e h[n]

    txSig = txSig + filterOut;    % sinal a ser transmitido     
end

% Registro do PAPR do sinal UFMC pr�-codificado
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);    % objeto pm para medi��o do PAPR
paprUFMC = pm(txSig);    % registro do PAPR do sinal UFMC
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);    % mostrar na tela

% Registro do PAPR do sinal OFDM pr�-codificado
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);     % objeto pm para medi��o do PAPR
paprOFDM = pm(ifftOut);    % registro do PAPR do sinal OFDM
disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);    % mostrar na tela