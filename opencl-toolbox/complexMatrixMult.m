clear all; clc;

% Inicialização e compilação do kernel
ocl = opencl();
ocl.initialize();
ocl.addfile('cl/complex_matrix_multiply.cl');
ocl.build();

% Definição das matrizes complexas
matrixA = int32([1+2i, 3+4i; 5+6i, 7+8i]); % declaração da matriz A
M = size(matrixA, 1); % linhas da matriz A
K = size(matrixA, 2); % colunas da matriz A

matrixB = int32([9+10i, 11+12i; 13+14i, 15+16i]); % declaração da matriz B
N = size(matrixB, 2); % linhas da matriz B

matrixC = zeros(M, N, 'int32'); % matriz C[M, N]

% Dividindo as partes reais e imaginárias das matrizes A e B
A_real = real(matrixA);
A_imag = imag(matrixA);
B_real = real(matrixB);
B_imag = imag(matrixB);

r_matrixA = clobject(reshape(A_real', 1, []));
i_matrixA = clobject(reshape(A_imag', 1, []));
r_matrixB = clobject(reshape(B_real', 1, []));
i_matrixB = clobject(reshape(B_imag', 1, []));
r_matrixC = clobject(matrixC);
i_matrixC = clobject(matrixC);

% Configuração dos tamanhos globais e locais do trabalho
global_work_size = [N, M, 0];
local_work_size = [1, 1, 0];

% Compilação e execução do kernel
kernel = clkernel('complex_matrix_multiply', global_work_size, local_work_size);
kernel(r_matrixA, i_matrixA, r_matrixB, i_matrixB, r_matrixC, i_matrixC, int32(M), int32(N), int32(K));

ocl.wait();

r_A = reshape(r_matrixA.get(), N, M)';
i_A = reshape(i_matrixA.get(), N, M)';
r_B = reshape(r_matrixB.get(), N, M)';
i_B = reshape(i_matrixB.get(), N, M)';
r_C = reshape(r_matrixC.get(), N, M)';
i_C = reshape(i_matrixC.get(), N, M)';

A = complex(r_A, i_A);
B = complex(r_B, i_B);
C = complex(r_C, i_C);

disp(A);
disp(B);
disp(C);
