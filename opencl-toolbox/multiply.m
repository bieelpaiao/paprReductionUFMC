% Limpar workspace e console
clear all
clc

% Inicialização e compilação do kernel OpenCL
ocl = opencl();
ocl.initialize();
ocl.addfile('cl/multiply.cl'); % Adicionar o arquivo do kernel
ocl.build();

% Definir os números inteiros a serem multiplicados
a = 8;
b = 9;

% Criar objetos clobject para os números inteiros
a_buffer = clobject(int32(a));
b_buffer = clobject(int32(b));

% Definir array para armazenar o resultado da multiplicação
result_buffer = clobject(zeros(1, 1, 'int32'));

% Definir as dimensões globais e locais do trabalho
global_work_size = [1, 0, 0]; % Apenas um work-item
local_work_size = [1, 0, 0]; % Um work-item por grupo local

% Criar objeto clkernel para o kernel de multiplicação de inteiros
multiply_kernel = clkernel('multiplyIntegers', global_work_size, local_work_size);

% Executar o kernel passando os números inteiros como argumentos
multiply_kernel(a_buffer, b_buffer, result_buffer);

% Esperar pela conclusão da execução do kernel
ocl.wait();

% Recuperar o resultado da multiplicação
result = result_buffer.get();

% Exibir o resultado
disp(['Resultado da multiplicação: ', num2str(result)]);
