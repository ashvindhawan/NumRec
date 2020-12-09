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

#define UNROLL 8
// #define BLOCKSIZE (64 / sizeof (double))
#define BLOCKSIZE 64

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
        // double * this_row = malloc(cols*sizeof(double));
        double * this_row = calloc(cols, sizeof(double));
        // for(int col = 0; col<cols; col++) {
        //     *(this_row+col) = 0;
        // }
        *(data+row) = this_row;
    }
    (*mat)->data = data;

    // double * act_data = calloc(cols*rows, sizeof(double));

    // #pragma omp parallel for
    // for(int row = 0; row < rows; row++) {
    //     data[row] = &(act_data[row * cols]);
    // }
    // (*mat)->data = data;
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

    (*mat)->ref_cnt = 1;
    matrix* curParent = (*mat)->parent;
    while(curParent!=NULL) {
        curParent->ref_cnt = curParent->ref_cnt + 1;
        curParent = curParent->parent;
    }
    //from->ref_cnt = from->ref_cnt + 1; this needs to bubble all the way up 
    
    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }
    double ** data = malloc(rows*sizeof(double*));

    // #pragma omp parallel for
    for (int i = 0; i < rows ; i++) {       
        double *addr = from->data[row_offset+i] + col_offset;  
        data[i] = addr;    
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
    // if(mat == NULL) {
    //     return;
    // } else if (mat->ref_cnt <= 1) {
    //     int rows = mat->rows;
    //     for(int row = 0; row<rows; row++) {
    //         free((mat->data)[row]);
    //     }
    //     free(mat);
    //     return;
    // } else {
    //     free(mat);
    // }
    if (mat == NULL) {
        return;
    }
    mat->ref_cnt = mat->ref_cnt - 1;
    if(mat->parent != NULL) {
        deallocate_matrix(mat->parent);
    }
    if(mat->ref_cnt==0) {
	    if(mat->parent==NULL) {
        	for(int row = 0; row<mat->rows; row++) {
                double ** data = mat->data;
            	free(data[row]);
            }  
            // free(*(mat->data));
            free(mat->data);
        }
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
    #pragma omp parallel for collapse(2)
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

    // #pragma omp parallel for collapse(2)
    // for (int i = 0; i<rows; i++) {
    //     for (int j = 0; j<cols; j++) {
    //         double elem1 = get(mat1, i, j);
    //         double elem2 = get(mat2, i, j);
    //         set(result, i, j, elem1 + elem2);
    //     }
    // }

    double **A = mat1->data;
    double **B = mat2->data;
    double **C = result->data;

    #pragma omp parallel for
    for (int i = 0; i<rows; ++i) {
        memcpy(C[i], B[i], cols * sizeof(double));
        for (int j = 0; j<cols / 4 * 4; j+= 4) {
            // _mm256_storeu_pd(C[i] + j, _mm256_add_pd (_mm256_loadu_pd (A[i] + j), _mm256_loadu_pd (B[i]+ j)));
            _mm256_storeu_pd(C[i] + j, _mm256_add_pd (_mm256_loadu_pd (A[i] + j), _mm256_loadu_pd (C[i] + j)));
        }

        for (int j = cols / 4 * 4; j < cols; ++j) {
            C[i][j] = C[i][j] + B[i][j];
        }
    }
    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rows = result->rows;
    int cols = result->cols;
    // for (int i = 0; i<rows; i++) {
    //     for (int j = 0; j<cols; j++) {
    //         double elem1 = get(mat1, i, j);
    //         double elem2 = get(mat2, i, j);
    //         set(result, i, j, elem1 - elem2);
    //     }
    // }

    double **A = mat1->data;
    double **B = mat2->data;
    double **C = result->data;

    #pragma omp parallel for
    for (int i = 0; i<rows; ++i) {
        for (int j = 0; j<cols / 4 * 4; j+= 4) {
            _mm256_storeu_pd(C[i] + j, _mm256_sub_pd (_mm256_loadu_pd (A[i] + j), _mm256_loadu_pd (B[i]+ j)));
        }

        for (int j = cols / 4 * 4; j < cols; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return 0;
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

#define UNROLL 8
// #define BLOCKSIZE (64 / sizeof (double))
#define BLOCKSIZE 64


int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) { 
    //Henry: need to add remainders (i.e. when matrix isn't a multiple of blocksize) otherwise it segfaults
    /* TODO: YOUR CODE HERE */
    // int rows = result->rows;
    // int cols = result->cols;
    // double sum;

    double ** mat1_data = mat1->data;
    double ** mat2_data = mat2->data;
    double ** result_data = result->data;
    
    int m1r = mat1->rows;
    int m1c = mat1->cols;
    int m2r = mat2->rows;
    int m2c = mat2->cols;

    // #pragma omp parallel for collapse(2)
    // for (int i = 0; i < m1c; i++) {
    //     for (int j = 0; j < m2r; j+=4*UNROLL) {
            
    //         __m256d c[UNROLL];

    //         for (int x = 0; x < UNROLL; x ++) {
    //             c[x] = _mm256_loadu_pd(result_data[i] + (j + x * 4));
    //         }
    //         for (int k = 0; k < m1r; k++) {
    //             __m256d a = _mm256_broadcast_sd(mat1_data[i] + k);
    //             for (int x = 0; x < UNROLL; x++) {
    //                 __m256d b = _mm256_loadu_pd(mat2_data[k]+ (j+4*x));
    //                 c[x] = _mm256_fmadd_pd(a, b, c[x]);
    //             }
    //         }
            
    //         for (int x = 0; x < UNROLL; x++) {
    //             _mm256_storeu_pd(result_data[i]+ (j + x * 4), c[x]);
    //         }
    //     }
    // }             

    double * A = *mat1->data;
    double * B = *mat2->data;
    double * C = *result->data;
    // memset(C, 0, m1r * m2c *sizeof(double));
    #pragma omp parallel for collapse(2)
    for (int si = 0; si < m1c; si = si + BLOCKSIZE) {
        for (int sj = 0; sj < m2r; sj = sj + BLOCKSIZE) {
            // for (int i = si; i < si+BLOCKSIZE; ++i) {
            //     for (int j = sj; j < sj + BLOCKSIZE; ++j) {
            //         memset(&result_data[i][j], 0, BLOCKSIZE *sizeof(double));
            //     }
            // }
            memset(&result_data[si][sj], 0, BLOCKSIZE *sizeof(double)); // set result data to 0 
            for (int sk = 0; sk < m2c; sk = sk + BLOCKSIZE) { // iterate through matrix 2's columns
                for (int i=si; i < si+BLOCKSIZE; ++i) { //loop unrolling
                    _mm_prefetch (&result_data[i][sj + BLOCKSIZE / 2], _MM_HINT_NTA);
                    // _mm_prefetch (C + i * m2c + (sj + BLOCKSIZE / 2), _MM_HINT_NTA);
                    for (int j =sj; j < sj+BLOCKSIZE; j = j + 4*UNROLL) {
                        __m256d c[UNROLL];
                        for (int x = 0; x < UNROLL; ++x) {
                            c[x] = _mm256_loadu_pd(&result_data[i][j + x * 4]); //creating a place to store cumulative sum
                            // c[x] = _mm256_loadu_pd(C + i * m2c + (j + x * 4));
                        }
                        for (int k = sk; k < sk+BLOCKSIZE; ++k) {
                            __m256d a = _mm256_broadcast_sd(&mat1_data[i][k]);
                            // __m256d a = _mm256_broadcast_sd(A + i * m1c + k);
                            for (int x = 0; x < UNROLL; ++x) {
                                __m256d b = _mm256_loadu_pd(&mat2_data[k][j+4*x]);
                                // __m256d b = _mm256_loadu_pd(B + k * m2c + (j+4*x));
                                c[x] = _mm256_fmadd_pd(a, b, c[x]);
                            }
                        }
                        for (int x = 0; x < UNROLL; ++x) {
                            _mm256_storeu_pd(&result_data[i][j + x * 4], c[x]);
                            // _mm256_storeu_pd(C + i * m2c + (j + x * 4), c[x]);
                        }
                    }
                }
            }
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
    // mul_matrix(*placeholder, mat, mat);
    #pragma omp parallel for
    for(int i = 0; i < mat->rows; i++) {
        set(result, i, i, 1.0);
        set(*placeholder, i, i, 1.0);
    }

    // for (int i = 0; i < pow; i++) {
    //     mul_matrix(result, *placeholder, mat);
    //     copy(result, *placeholder);
    // }

    matrix * t1 = result;
    matrix * t2 = *placeholder;

    // for (int i = 0; i < pow; i++) {
    //     mul_matrix(t1, t2, mat);
    //     matrix * temp = t1;
    //     t1 = t2;
    //     t2 = temp;
    // }

    // if (pow % 2 == 0) {
    //     copy(*placeholder, result);
    // }
    
    // deallocate_matrix(*placeholder);
    // free(placeholder);
    // return 0;

    int swaps = 0;
    while (pow) {
        if (pow % 2 == 1) {
            mul_matrix(t1, t2, mat);
            pow -= 1;
            matrix * temp = t1;
            t1 = t2;
            t2 = temp;
        } else {
            mul_matrix(t1, t2, t2);
            pow = pow / 2;
            matrix * temp = t1;
            t1 = t2;
            t2 = temp;
        }
        ++swaps;
    }

    if (swaps % 2 == 0) {
        copy(*placeholder, result);
    }

    deallocate_matrix(*placeholder);
    free(placeholder);
    return 0;
}


void copy(matrix * src, matrix* dest) {
    // int** src_data = src->data;
    // int** dest_data = dest->data;
    int rows = src->rows;
    int cols = src->cols;
    #pragma omp parallel for collapse(2)
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

