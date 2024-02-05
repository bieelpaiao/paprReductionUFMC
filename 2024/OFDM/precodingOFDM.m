% clear all; clc;

%%
M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 256;           % Number of OFDM subcarriers
cpLen = 32;            % OFDM cyclic prefix length
nSym = 100;
scs = 1000000;

s = rng(211);          % configuração do gerador de números aleatórios

%%
ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen, 'NumGuardBandCarriers', [6;6], 'NumSymbols', nSym);
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen, 'NumGuardBandCarriers', [6;6], 'NumSymbols', nSym);

%%
channel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

%%
ofdmDims = info(ofdmMod);

%%
beta = 0.35;
inputSize = ofdmDims.DataInputSize(1);

%%
numDC = round(inputSize/(1+beta));
Np = inputSize-numDC;

%%
frameSize = [k*(inputSize-Np) nSym];

%%
pim = matrixGeneration(frameSize(1)/k, Np);

%%
EbNo = 10;
snr = EbNo + 10*log10(k) + 10*log10(numDC/numSC);

%%
dataIn = randi([0,1],frameSize);                                % Generate binary data
qamTx = qammod(dataIn, M, "gray", "InputType", "bit");          % Apply QPSK modulation
qamFiltered = precod(pim, qamTx);                               % Apply SRC precoding
txSig = ofdmMod(qamFiltered);                                   % Apply OFDM modulation
powerDB = 10*log10(var(txSig));                                 % Calculate Tx signal power
noiseVar = 10.^(0.1*(powerDB-snr));                             % Calculate the noise variance
rxSig = channel(txSig,noiseVar);                                % Pass the signal through a noisy channel
qamRx = ofdmDemod(rxSig);                                       % Apply OFDM demodulation
qamDecoded = decod(pim, qamRx);                                 % Apply SRC* precoding
dataOut = qamdemod(qamDecoded, M, "gray", "OutputType","bit");  % Apply QPSK demodulation

%%
% pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
% paprOFDM = pm(txSig);
% disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);
% plotCCDF(pm);

%% Plot Spectrum Analyzer
Fs = ofdmMod.FFTLength * scs * ofdmMod.OversamplingFactor;
spectrum = spectrumAnalyzer('SampleRate', Fs,'ChannelNames', {'Tradicional', 'Pré-codificado'}, 'Title', 'Espectro OFDM');
% spectrum(txSig);




