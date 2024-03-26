// Kernel para multiplicação de dois números inteiros

__kernel void multiplyIntegers(__global int* a,
                                __global int* b,
                                __global int* result) {
    // Obtenção do índice global do work-item
    int index = get_global_id(0);

    // Multiplicação dos números e armazenamento no array de resultados
    result[index] = a[index] * b[index];
}
