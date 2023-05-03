#include "compute.h"

// Computes the dot product of vec1 and vec2, both of size n
int dot(uint32_t n, int32_t *vec1, int32_t *vec2) {
  // TODO: implement dot product of vec1 and vec2, both of size n
  int dot=0; 
  for(int i=0; i<n;i++){
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

  for (int32_t i=0; i< b_len; i++) {
	  vec_b_flipped[i] = b_matrix->data[b_len-i-1];
  }

  

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
