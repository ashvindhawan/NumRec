#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    if (rows<=0 || cols <=0) {
        PyErr_SetString(PyExc_ValueError, "Incorrect values for row and col");
        return -1;
    }
    *mat = (matrix *) malloc(sizeof(matrix));
    if (*mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "malloc failed!");
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->parent = NULL;
    (*mat)->ref_cnt = 1;
    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }
    double ** data = malloc(rows*sizeof(double*));
    for(int row = 0; row<rows; row++) {
        double * this_row = malloc(cols*sizeof(double));
        for(int col = 0; col<cols; col++) {
            *(this_row+col) = 0;
        }
        *(data+row) = this_row;
    }
    (*mat)->data = data;
    return 0;

}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows<=0 || cols <=0) {
        PyErr_SetString(PyExc_ValueError, "Incorrect values for row and col");
        return -1;
    }
    *mat = (matrix *) malloc(sizeof(matrix));
    if (*mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "malloc failed!");
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->parent = from;
    
    from->ref_cnt = from->ref_cnt + 1;
    (*mat)->ref_cnt = 2;
    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }
    double ** data = malloc(rows*sizeof(double*));
    for(int row = 0; row<rows; row++) {
        data[row] = &from->data[row+row_offset][col_offset];
    }
    (*mat)->data = data;
    return 0;
    
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    if(mat == NULL) {
        return;
    } else if (mat->ref_cnt <= 1) {
        int rows = mat->rows;
        for(int row = 0; row<rows; row++) {
            free((mat->data)[row]);
        }
        free(mat);
        return;
    } else {
        free(mat);
    }
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    double **data = mat->data;
    return (double) data[row][col];
    /* TODO: YOUR CODE HERE */
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    double** data = mat->data; 
    data[row][col] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    int row  = mat->rows;
    int col = mat->cols;
    for (int i = 0; i<row; i++) {
        for (int j = 0; j<col; j++) {
            set(mat, i, j, val);
        }
    }
    
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rows = result->rows;
    int cols = result->cols;
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<cols; j++) {
            double elem1 = get(mat1, i, j);
            double elem2 = get(mat2, i, j);
            set(result, i, j, elem1 + elem2);
        }
    }
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rows = result->rows;
    int cols = result->cols;
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<cols; j++) {
            double elem1 = get(mat1, i, j);
            double elem2 = get(mat2, i, j);
            set(result, i, j, elem1 - elem2);
        }
    }
}

double dot(double* row, double* col, int len) {
    double ans = 0;
    for(int i = 0; i < len; i ++) {
        ans = ans + (row[i] * col[i]);
    }
    return ans;
}

double* column(matrix *mat, int index) {
    double * col = malloc((mat->rows)*sizeof(int));
    for(int i = 0; i < mat->rows; i++) {
        col[i] = get(mat, i, index);
    }
    return col;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rows = result->rows;
    int cols = result->cols;
    double ** mat1_data = mat1->data;
    double ** mat2_data = mat2->data;
    for (int row = 0; row<rows; row++) {
        for(int col = 0; col<cols; col++) {
            double val = dot(mat1_data[row], column(mat2, col), cols);
            set(result, row, col, val);
        }
    }
    return 0;
}


/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    matrix ** placeholder = malloc(sizeof(matrix *));
    allocate_matrix(placeholder, mat->rows, mat->cols);
    mul_matrix(*placeholder, mat, mat);
    for (int i = 0; i < pow-2; i++) {
        mul_matrix(result, *placeholder, mat);
        copy(result, *placeholder);
    }
    return 0;
}

void copy(matrix * src, matrix* dest) {
    // int** src_data = src->data;
    // int** dest_data = dest->data;
    int rows = src->rows;
    int cols = src->cols;
    for(int row = 0; row < rows; row++) {
        for(int col = 0; col < cols; col++) {
            set(dest, row, col, get(src, row, col));
        }
    }
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    int rows = result->rows;
    int cols = result->cols;
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<cols; j++) {
            double elem = get(mat, i, j);
            set(result, i, j, -elem);
        }
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    int rows = result->rows;
    int cols = result->cols;
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<cols; j++) {
            double elem = get(mat, i, j);
            if (elem < 0) {
                elem = elem * -1;
            }
            set(result, i, j, elem);
        }
    }
    return 0;
}

