__kernel void complex_matrix_mult(__global double* A_real,
                                  __global double* A_imag,
                                  __global double* B_real,
                                  __global double* B_imag,
                                  __global double* C_real,
                                  __global double* C_imag,
                                  const int M, 
                                  const int N, 
                                  const int K) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    // Inicializa os valores da parte real e imagin�ria do elemento resultante como 0
    C_real[i * K + j] = 0;
    C_imag[i * K + j] = 0;

    // Realiza a multiplica��o de cada elemento da linha i de A pela coluna j de B
    for (int k = 0; k < N; k++) {
        // Calcula a parte real e imagin�ria do produto
        C_real[i * K + j] += A_real[i * N + k] * B_real[k * K + j] - A_imag[i * N + k] * B_imag[k * K + j];
        C_imag[i * K + j] += A_real[i * N + k] * B_imag[k * K + j] + A_imag[i * N + k] * B_real[k * K + j];
    }
}
