// Kernel para multiplicação de matrizes com números complexos
__kernel void complex_matrix_multiply(__global float2* A_real,
                                      __global float2* A_imag,
                                      __global float2* B_real,
                                      __global float2* B_imag,
                                      __global float2* C_real,
                                      __global float2* C_imag,
                                      int M, int N, int K) {
    int row = get_global_id(0);
    int col = get_global_id(1);

    float2 sum = (float2)(0.0f, 0.0f);

    for (int k = 0; k < K; ++k) {
        float2 a_real = A_real[row * K + k];
        float2 a_imag = A_imag[row * K + k];
        float2 b_real = B_real[k * N + col];
        float2 b_imag = B_imag[k * N + col];

        // Multiplicação complexa
        sum += (float2)(a_real.x * b_real.x - a_imag.x * b_imag.x - a_real.y * b_imag.y - a_imag.y * b_real.y,
                        a_real.x * b_imag.x + a_imag.x * b_real.x + a_real.y * b_real.y - a_imag.y * b_imag.y);
    }

    C_real[row * N + col] = sum;
    C_imag[row * N + col] = (float2)(0.0f, 0.0f); // Parte imaginária é sempre zero na multiplicação real
}
