clear all; clc;

% Inicializar o OpenCL
ocl = opencl();
ocl.initialize();

% Criar um buffer para a matriz 3x2
matriz = clbuffer('rw', 'uint64', 3*2);

% Definir a matriz 3x2
matriz_valores = [1 2; 3 4; 5 6];

% Escrever a matriz no buffer OpenCL
matriz.set(matriz_valores);

values = matriz.get();

% Limpar o buffer quando não for mais necessário
clear matriz;
