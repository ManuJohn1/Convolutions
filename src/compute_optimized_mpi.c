#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the dot product of vec1 and vec2, both of size n
int dot(uint32_t n, int32_t *vec1, int32_t *vec2) {
    // TODO: implement dot product of vec1 and vec2, both of size n
    int32_t sum = 0;
    for(uint32_t i = 0; i < n/64 * 64; i += 64) {
      __m256i vec_1a = _mm256_loadu_si256((__m256i *) (vec1 + i));
      __m256i vec_2a = _mm256_loadu_si256((__m256i *) (vec2 + i));
      __m256i sum_1 = _mm256_mullo_epi32(vec_1a, vec_2a);


      int32_t tmp_arr[8];
      _mm256_storeu_si256((__m256i *) tmp_arr, sum_1);
      sum += tmp_arr[0] + tmp_arr[1] + tmp_arr[2] + tmp_arr[3] + tmp_arr[4] + tmp_arr[5] + tmp_arr[6] + tmp_arr[7];

      __m256i vec_1b = _mm256_loadu_si256((__m256i *) (vec1 + i + 8));
      __m256i vec_2b = _mm256_loadu_si256((__m256i *) (vec2 + i + 8));
      __m256i sum_2 = _mm256_mullo_epi32(vec_1b, vec_2b);

      int32_t tmp_arr2[8];
      _mm256_storeu_si256((__m256i *) tmp_arr2, sum_2);
      sum += tmp_arr2[0] + tmp_arr2[1] + tmp_arr2[2] + tmp_arr2[3] + tmp_arr2[4] + tmp_arr2[5] + tmp_arr2[6] + tmp_arr2[7];

      __m256i vec_1c = _mm256_loadu_si256((__m256i *) (vec1 + i + 16));
      __m256i vec_2c = _mm256_loadu_si256((__m256i *) (vec2 + i + 16));
      __m256i sum_3 = _mm256_mullo_epi32(vec_1c, vec_2c);

      int32_t tmp_arr3[8];
      _mm256_storeu_si256((__m256i *) tmp_arr3, sum_3);
      sum += tmp_arr3[0] + tmp_arr3[1] + tmp_arr3[2] + tmp_arr3[3] + tmp_arr3[4] + tmp_arr3[5] + tmp_arr3[6] + tmp_arr3[7];

      __m256i vec_1d = _mm256_loadu_si256((__m256i *) (vec1 + i + 24));
      __m256i vec_2d = _mm256_loadu_si256((__m256i *) (vec2 + i + 24));
      __m256i sum_4 = _mm256_mullo_epi32(vec_1d, vec_2d);

      int32_t tmp_arr4[8];
      _mm256_storeu_si256((__m256i *) tmp_arr4, sum_4);
      sum += tmp_arr4[0] + tmp_arr4[1] + tmp_arr4[2] + tmp_arr4[3] + tmp_arr4[4] + tmp_arr4[5] + tmp_arr4[6] + tmp_arr4[7];

      __m256i vec_1e = _mm256_loadu_si256((__m256i *) (vec1 + i + 32));
      __m256i vec_2e = _mm256_loadu_si256((__m256i *) (vec2 + i + 32));
      __m256i sum_5 = _mm256_mullo_epi32(vec_1e, vec_2e);


      int32_t tmp_arr5[8];
      _mm256_storeu_si256((__m256i *) tmp_arr5, sum_5);
      sum += tmp_arr5[0] + tmp_arr5[1] + tmp_arr5[2] + tmp_arr5[3] + tmp_arr5[4] + tmp_arr5[5] + tmp_arr5[6] + tmp_arr5[7];

      __m256i vec_1f = _mm256_loadu_si256((__m256i *) (vec1 + i + 40));
      __m256i vec_2f = _mm256_loadu_si256((__m256i *) (vec2 + i + 40));
      __m256i sum_6 = _mm256_mullo_epi32(vec_1f, vec_2f);

      int32_t tmp_arr6[8];
      _mm256_storeu_si256((__m256i *) tmp_arr6, sum_6);
      sum += tmp_arr6[0] + tmp_arr6[1] + tmp_arr6[2] + tmp_arr6[3] + tmp_arr6[4] + tmp_arr6[5] + tmp_arr6[6] + tmp_arr6[7];

      __m256i vec_1g = _mm256_loadu_si256((__m256i *) (vec1 + i + 48));
      __m256i vec_2g = _mm256_loadu_si256((__m256i *) (vec2 + i + 48));
      __m256i sum_7 = _mm256_mullo_epi32(vec_1g, vec_2g);

      int32_t tmp_arr7[8];
      _mm256_storeu_si256((__m256i *) tmp_arr7, sum_7);
      sum += tmp_arr7[0] + tmp_arr7[1] + tmp_arr7[2] + tmp_arr7[3] + tmp_arr7[4] + tmp_arr7[5] + tmp_arr7[6] + tmp_arr7[7];
        
      __m256i vec_1h = _mm256_loadu_si256((__m256i *) (vec1 + i + 56));
      __m256i vec_2h = _mm256_loadu_si256((__m256i *) (vec2 + i + 56));
      __m256i sum_8 = _mm256_mullo_epi32(vec_1h, vec_2h);

      int32_t tmp_arr8[8];
      _mm256_storeu_si256((__m256i *) tmp_arr8, sum_8);
      sum += tmp_arr8[0] + tmp_arr8[1] + tmp_arr8[2] + tmp_arr8[3] + tmp_arr8[4] + tmp_arr8[5] + tmp_arr8[6] + tmp_arr8[7];
}

    int temp = n/64*64;
    for(unsigned int i = temp; i < n/8*8; i+=8) {
      __m256i vec_1a = _mm256_loadu_si256((__m256i *) (vec1 + i));
      __m256i vec_2a = _mm256_loadu_si256((__m256i *) (vec2 + i));
      __m256i sum_1 = _mm256_mullo_epi32(vec_1a, vec_2a);


      int32_t tmp_arr[8];
      _mm256_storeu_si256((__m256i *) tmp_arr, sum_1);
      sum += tmp_arr[0] + tmp_arr[1] + tmp_arr[2] + tmp_arr[3] + tmp_arr[4] + tmp_arr[5] + tmp_arr[6] + tmp_arr[7];

    }

    for(unsigned int i = n/8*8; i < n; i++) {
        sum += vec1[i]*vec2[i];
    }
  return sum;
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
