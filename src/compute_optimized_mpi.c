#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the dot product of vec1 and vec2, both of size n
int dot(uint32_t n, int32_t *vec1, int32_t *vec2) {
  // TODO: implement dot product of vec1 and vec2, both of size n
  int32_t dot = 0;
  for(uint32_t i = 0; i < n/16 * 16; i += 16) {
      __m256i vec_1a = _mm256_loadu_si256((__m256i *) (vec1 + i));
      __m256i vec_2a = _mm256_loadu_si256((__m256i *) (vec2 + i));
      __m256i sum_1 = _mm256_mullo_epi32(vec_1a, vec_2a);


      int32_t tmp_arr[8];
      _mm256_storeu_si256((__m256i *) tmp_arr, sum_1);
      dot += tmp_arr[0] + tmp_arr[1] + tmp_arr[2] + tmp_arr[3] + tmp_arr[4] + tmp_arr[5] + tmp_arr[6] + tmp_arr[7];

      __m256i vec_1b = _mm256_loadu_si256((__m256i *) (vec1 + i + 8));
      __m256i vec_2b = _mm256_loadu_si256((__m256i *) (vec2 + i + 8));
      __m256i sum_2 = _mm256_mullo_epi32(vec_1b, vec_2b);

      int32_t tmp_arr2[8];
      _mm256_storeu_si256((__m256i *) tmp_arr2, sum_2);
      dot += tmp_arr2[0] + tmp_arr2[1] + tmp_arr2[2] + tmp_arr2[3] + tmp_arr2[4] + tmp_arr2[5] + tmp_arr2[6] + tmp_arr2[7];
    }

    for(unsigned int i = n/16 * 16; i < n; i++) {
        dot += vec1[i]*vec2[i];
    }
  return dot;
}

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  // output_matrix
  matrix_t *res_matrix = malloc(sizeof(matrix_t));
  res_matrix->rows = a_matrix->rows - b_matrix->rows + 1;
  res_matrix->cols = a_matrix->cols - b_matrix->cols + 1;
  res_matrix->data = malloc(sizeof(int32_t) * res_matrix->rows*res_matrix->cols);


  int32_t b_len = b_matrix->rows*b_matrix->cols;
  int32_t *vec_b_flipped = malloc(sizeof(int32_t)*b_len);
  
  #pragma omp for
  for (int32_t i=0; i< b_len; i++) {
      vec_b_flipped[i] = b_matrix->data[b_len-i-1];
  }
  

  #pragma omp parallel for collapse(2)
      for (int32_t p=0; p < res_matrix->rows; p++) {
          for (int32_t q=0; q < res_matrix->cols; q++) {
              int32_t sum = 0;
              for(int32_t i=0; i < b_matrix->rows; i++) {
                  sum += dot(b_matrix->cols, (a_matrix->data + (p+i)*a_matrix->cols + q), (vec_b_flipped + (i)*b_matrix->cols));
              }
              res_matrix->data[p*res_matrix->cols + q] = sum;
          } 
      }
  free(vec_b_flipped);
  *output_matrix = res_matrix;
  return 0;
}

// Executes a task
int execute_task(task_t *task) {
  matrix_t *a_matrix, *b_matrix, *output_matrix;

  if (read_matrix(get_a_matrix_path(task), &a_matrix))
    return -1;
  if (read_matrix(get_b_matrix_path(task), &b_matrix))
    return -1;

  if (convolve(a_matrix, b_matrix, &output_matrix))
    return -1;

  if (write_matrix(get_output_matrix_path(task), output_matrix))
    return -1;

  free(a_matrix->data);
  free(b_matrix->data);
  free(output_matrix->data);
  free(a_matrix);
  free(b_matrix);
  free(output_matrix);
  return 0;
}
