% Parâmetros
fs = 10000;                    % Frequência de amostragem (Hz) para boa resolução
t = 0.097:1/fs:0.108;          % Vetor de tempo ajustado para 5 ciclos da maior frequência próximo ao ponto de fase

% Frequências das cinco senóides
frequencies = [50, 120, 200, 300, 400];  

% Inicializa o sinal resultante
resultant_signal = zeros(size(t));

% Cria uma nova figura
figure;
hold on;

% Gera cada senóide e as plota
for i = 1:length(frequencies)
    sine_wave = cos(2 * pi * frequencies(i) * t);
    resultant_signal = resultant_signal + sine_wave;  % Soma cada senóide ao sinal resultante
    plot(t, sine_wave, 'DisplayName', ['Senoide ' num2str(i) ' (' num2str(frequencies(i)) ' Hz)'], 'LineWidth', 2);
end

% Plota o sinal resultante (soma de todas as senóides)
plot(t, resultant_signal, 'k', 'LineWidth', 2, 'DisplayName', 'Soma das Senoides');

% Adiciona uma linha vertical pontilhada em t = 20 ms
xline(0.1, 'r--');

% Configurações do gráfico
title('Demonstração de Alta PAPR pela Soma de Senoides');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;
hold off;
ylim([-3 6]);
xlim([0.097 0.108]);
set(gca, 'LineWidth', 1.5);
