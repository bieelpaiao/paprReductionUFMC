clear all: clc

%%
M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 512;           % Number of OFDM subcarriers
cpLen = 32;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted

% s = rng(211);       % configuração do gerador de números aleatórios

%%
ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);

%%
channel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

%%
errorRate = comm.ErrorRate('ResetInputPort',true);

%%
ofdmDims = info(ofdmMod);

%%
numDC = ofdmDims.DataInputSize(1);

%%
frameSize = [k*numDC 1];

%%
EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC);

%%
berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);

%%
for m = 1:length(EbNoVec)
    snr = snrVec(m);
    
    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        dataIn = randi([0,1],frameSize);                            % Generate binary data
        qamTx = qammod(dataIn, M, "gray", "InputType", "bit");      % Apply QPSK modulation
        txSig = ofdmMod(qamTx);                                     % Apply OFDM modulation
        powerDB = 10*log10(var(txSig));                             % Calculate Tx signal power
        noiseVar = 10.^(0.1*(powerDB-snr));                         % Calculate the noise variance
        rxSig = channel(txSig,noiseVar);                            % Pass the signal through a noisy channel
        qamRx = ofdmDemod(rxSig);                                   % Apply OFDM demodulation
        dataOut = qamdemod(qamRx, M, "gray", "OutputType","bit");   % Apply QPSK demodulation
        errorStats = errorRate(dataIn,dataOut,0);                   % Collect error statistics
    end
    
    berVec(m,:) = errorStats;                         % Save BER data
    errorStats = errorRate(dataIn,dataOut,1);         % Reset the error rate calculator
end

%%
berTheory = berawgn(EbNoVec,'qam',M,'nondiff');

%%
figure (1)
% semilogy(EbNoVec,berVec(:,1),'*')
semilogy(EbNoVec,berVec(:,1), ':o')
% hold on
% semilogy(EbNoVec,berTheory)
% legend('Simulation','Theory','Location','Best')
% legend('Simulation')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off
