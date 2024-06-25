% clear all; clc;

% Inicialização e compilação do kernel
ocl = opencl();
ocl.initialize();
ocl.addfile('cl/matrix_multiplication.cl');
ocl.build();

matrixA = int32(reshape(1:10000, 100, 100)); % declaração da matriz A
M = size(matrixA, 1); % linhas da matriz A
K = size(matrixA, 2); % colunas da matriz A

matrixB = int32(matrixA'); % declaração da matriz B
N = size(matrixB, 2); % linhas da matriz B

matrixC = zeros(M, N, 'int32'); % matriz C[M, N]

% Criação das matrizes de entrada
matrix1 = clobject(reshape(matrixA', 1, []));
matrix2 = clobject(reshape(matrixB', 1, []));
matrix3 = clobject(matrixC);

% Definição das dimensões globais e locais do trabalho
global_work_size = [N, M, 0];
local_work_size = [1, 1, 0];

% Criação do objeto clkernel para o kernel setOnes
setOnesKernel = clkernel('multiplyMatrices', global_work_size, local_work_size);

% Execução do kernel passando as matrizes de entrada como argumentos
setOnesKernel(matrix1, matrix2, matrix3, int32(M), int32(N), int32(K));

% Espera pela conclusão da execução do kernel
ocl.wait();

% Recuperação dos valores das matrizes de saída
resultMatrix1 = reshape(matrix1.get(), K, M)';
resultMatrix2 = reshape(matrix2.get(), N, K)';
resultMatrix3 = reshape(matrix3.get(), N, M)';

% Exibir os resultados
% disp('Matrix 1:');
% disp(resultMatrix1);
% disp('Matrix 2:');
% disp(resultMatrix2);
% disp('Matrix 3:');
% disp(resultMatrix3);
% disp(time_kernel)

delete(matrix1);
delete(matrix2);
delete(matrix3);
