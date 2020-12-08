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

    // for(int row = 0; row<rows; row++) {
    //     double * this_row = calloc(cols, sizeof(double));
    //     *(data+row) = this_row;
    // }
    // (*mat)->data = data;

    double * act_data = calloc(cols*rows, sizeof(double));

    #pragma omp parallel for
    for(int row = 0; row < rows; row++) {
        data[row] = &(act_data[row * cols]);
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
    if (mat == NULL) {
        return;
    }
    mat->ref_cnt = mat->ref_cnt - 1;
    if(mat->parent != NULL) {
        deallocate_matrix(mat->parent);
    }
    if(mat->ref_cnt==0) {
	    if(mat->parent==NULL) {
        	// for(int row = 0; row<mat->rows; row++) {
            //     double ** data = mat->data;
            // 	free(data[row]);
            // }  
            free(*(mat->data));
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
    return (*(mat->data))[row * mat->cols + col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    (*(mat->data))[row * mat->cols + col] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    int rows  = mat->rows;
    int cols = mat->cols;

    int boundary = cols / 4 * 4;
    double * data = *(mat->data);

    __m256d vals = _mm256_set1_pd(val);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<boundary; j+= 4) {
            int offset = i * cols + j;
            _mm256_storeu_pd(data + offset, vals);
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = boundary; j<cols; j++) {
            int offset = i * cols + j;
            data[offset] = val;
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

    double *A = *(mat1->data);
    double *B = *(mat2->data);
    double *C = *(result->data);

    int boundary = cols / 4 * 4;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<boundary; j+= 4) {
            int offset = i * cols + j;
            __m256d a = _mm256_loadu_pd(A + offset);
            __m256d b = _mm256_loadu_pd(B + offset);
            _mm256_storeu_pd(C + offset, _mm256_add_pd(a, b));
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = boundary; j<cols; j++) {
            int offset = i * cols + j;
            C[offset] = A[offset] + B[offset];
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

    double *A = *(mat1->data);
    double *B = *(mat2->data);
    double *C = *(result->data);

    int boundary = cols / 4 * 4;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<boundary; j+= 4) {
            int offset = i * cols + j;
            __m256d a = _mm256_loadu_pd(A + offset);
            __m256d b = _mm256_loadu_pd(B + offset);
            _mm256_storeu_pd(C + offset, _mm256_sub_pd(a, b));
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = boundary; j<cols; j++) {
            int offset = i * cols + j;
            C[offset] = A[offset] - B[offset];
        }
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */

#define UNROLL 4
#define BLOCKSIZE 64

int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    double ** mat1_data = mat1->data;
    double ** mat2_data = mat2->data;
    double ** result_data = result->data;
    
    int m1r = mat1->rows, m1c = mat1->cols, m2r = mat2->rows, m2c = mat2->cols;

    // int m2c_e = m2c / (4*UNROLL) * 4*UNROLL;

    // double * A = *(mat1->data);
    // double * B = *(mat2->data);
    // double * C = *(result->data);

    // #pragma omp parallel for collapse(2)
    // for (int i = 0; i < m1r; i++) {
    //     for (int j = 0; j < m2c_e; j+=4*UNROLL) {
    //         __m256d c[UNROLL];

    //         for (int x = 0; x < UNROLL; x ++)
    //             c[x] = _mm256_setzero_pd();

    //         for (int k = 0; k < m1c; k++) {
    //             __m256d a = _mm256_broadcast_sd(A + i * m1c + k);
    //             for (int x = 0; x < UNROLL; x++) {
    //                 __m256d b = _mm256_loadu_pd(B + k * m2c + (j+4*x));
    //                 c[x] = _mm256_fmadd_pd(a, b, c[x]);
    //             }
    //         }
            
    //         for (int x = 0; x < UNROLL; x++)
    //             _mm256_storeu_pd(C + i * m2c + (j + x * 4), c[x]);
    //     }
    // }    

    // #pragma omp parallel for collapse(2)
    // for (int i = 0; i < m1r; ++i) {
    //     for (int j = m2c_e; j < m2c; ++j) {
    //         double s = 0.0;
    //         for (int k = 0; k < m1c; ++k) {              
    //             s += A[i * m1c + k] * B[k * m2c + j];
    //         }
    //         C[i * m2c + j] = s;
    //     }
    // }         

    double * A = *mat1->data;
    double * B = *mat2->data;
    double * C = *result->data;

    int m1r_e = m1r / BLOCKSIZE * BLOCKSIZE;
    int m2r_e = m2r / BLOCKSIZE * BLOCKSIZE;
    int m2c_e = m2c / BLOCKSIZE * BLOCKSIZE;

    memset(*result_data, 0, m1r * m2c *sizeof(double));

    #pragma omp parallel for collapse(2)
    for (int si = 0; si < m1r_e; si = si + BLOCKSIZE) {
        for (int sj = 0; sj < m2c_e; sj = sj + BLOCKSIZE) { 
            for (int sk = 0; sk < m2r_e; sk = sk + BLOCKSIZE) {
                for (int i=si; i < si+BLOCKSIZE; ++i) {
                    _mm_prefetch (C + i * m2c + (sj + 0), _MM_HINT_NTA);
                    for (int j =sj; j < sj+BLOCKSIZE; j = j + 4*UNROLL) {
                        __m256d c[UNROLL];
                        for (int x = 0; x < UNROLL; ++x)
                            c[x] = _mm256_loadu_pd(C + i * m2c + j + x * 4);
                        for (int k = sk; k < sk+BLOCKSIZE; ++k) {
                            __m256d a = _mm256_broadcast_sd(A + i * m1c + k);
                            for (int x = 0; x < UNROLL; ++x) {
                                __m256d b = _mm256_loadu_pd(B + k * m2c + j+4*x);
                                c[x] = _mm256_fmadd_pd(a, b, c[x]);
                            }
                        }
                        for (int x = 0; x < UNROLL; ++x)
                            _mm256_storeu_pd(C + i * m2c + (j + x * 4), c[x]);
                    }
                }
            }

            for (int sk = m2r_e; sk < m2r; ++sk) {
                for (int i=si; i < si+BLOCKSIZE; ++i) {
                    _mm_prefetch (C + i * m2c + (sj + 0), _MM_HINT_NTA);
                    for (int j =sj; j < sj+BLOCKSIZE; j = j + 4*UNROLL) {
                        __m256d c[UNROLL];
                        for (int x = 0; x < UNROLL; ++x)
                            c[x] = _mm256_loadu_pd(C + i * m2c + j + x * 4);
                        
                        __m256d a = _mm256_broadcast_sd(A + i * m1c + sk);
                        for (int x = 0; x < UNROLL; ++x) {
                            __m256d b = _mm256_loadu_pd(B + sk * m2c + j+4*x);
                            c[x] = _mm256_fmadd_pd(a, b, c[x]);
                        }
                        for (int x = 0; x < UNROLL; ++x)
                            _mm256_storeu_pd(C + i * m2c + (j + x * 4), c[x]);
                    }
                }
            }
        }
    }

    for (int si = 0; si < m1r_e; si = si + BLOCKSIZE) {
        for (int sj = m2c_e; sj < m2c; sj = ++sj) { 
            for (int sk = 0; sk < m2r_e; sk = sk + BLOCKSIZE) {
                for (int i=si; i < si+BLOCKSIZE; ++i) {
                    double s = C[i * m2c + sj];
                    for (int k = sk; k < sk+BLOCKSIZE; ++k) {
                        s += A[i * m1c + k] * B[k * m2c + sj]; 
                    }
                    C[i * m2c + sj] = s;
                }
            }

            for (int sk = m2r_e; sk < m2r; ++sk) {
                for (int i=si; i < si+BLOCKSIZE; ++i) {
                    C[i * m2c + sj] += A[i * m1c + sk] * B[sk * m2c + sj];
                }
            }
        }
    }

    for (int si = m1r_e; si < m1r; si = ++si) {
        for (int sj = 0; sj < m2c_e; sj = sj + BLOCKSIZE) { 
            for (int sk = 0; sk < m2r_e; sk = sk + BLOCKSIZE) {
                _mm_prefetch (C + si * m2c + (sj + 0), _MM_HINT_NTA);
                for (int j =sj; j < sj+BLOCKSIZE; j = j + 4*UNROLL) {
                    __m256d c[UNROLL];
                    for (int x = 0; x < UNROLL; ++x)
                        c[x] = _mm256_loadu_pd(C + si * m2c + j + x * 4);
                    for (int k = sk; k < sk+BLOCKSIZE; ++k) {
                        __m256d a = _mm256_broadcast_sd(A + si * m1c + k);
                        for (int x = 0; x < UNROLL; ++x) {
                            __m256d b = _mm256_loadu_pd(B + k * m2c + j+4*x);
                            c[x] = _mm256_fmadd_pd(a, b, c[x]);
                        }
                    }
                    for (int x = 0; x < UNROLL; ++x)
                        _mm256_storeu_pd(C + si * m2c + (j + x * 4), c[x]);
                }
            }

            for (int sk = m2r_e; sk < m2r; ++sk) {
                _mm_prefetch (C + si * m2c + (sj + 0), _MM_HINT_NTA);
                for (int j =sj; j < sj+BLOCKSIZE; j = j + 4*UNROLL) {
                    __m256d c[UNROLL];
                    for (int x = 0; x < UNROLL; ++x)
                        c[x] = _mm256_loadu_pd(C + si * m2c + j + x * 4);
                    
                    __m256d a = _mm256_broadcast_sd(A + si * m1c + sk);
                    for (int x = 0; x < UNROLL; ++x) {
                        __m256d b = _mm256_loadu_pd(B + sk * m2c + j+4*x);
                        c[x] = _mm256_fmadd_pd(a, b, c[x]);
                    }
                    for (int x = 0; x < UNROLL; ++x)
                        _mm256_storeu_pd(C + si * m2c + (j + x * 4), c[x]);
                }
            }
        }
    }

    for (int si = m1r_e; si < m1r; si = ++si) {
        for (int sj = m2c_e; sj < m2c; sj = ++sj) { 
            for (int sk = 0; sk < m2r_e; sk = sk + BLOCKSIZE) {
                double s = C[si * m2c + sj];
                for (int k = sk; k < sk+BLOCKSIZE; ++k) {
                    s += A[si * m1c + k] * B[k * m2c + sj]; 
                }
                C[si * m2c + sj] = s;
            }

            for (int sk = m2r_e; sk < m2r; ++sk) {
                C[si * m2c + sj] += A[si * m1c + sk] * B[sk * m2c + sj];
            }
        }
    }
    return 0;
}

int mul_matrix_pow(matrix *result, matrix *mat1, matrix *mat2) {
    double ** mat1_data = mat1->data;
    double ** mat2_data = mat2->data;
    double ** result_data = result->data;
    
    int m1r = mat1->rows, m1c = mat1->cols, m2r = mat2->rows, m2c = mat2->cols;

    int m2c_e = m2c / (4*UNROLL) * 4*UNROLL;

    double * A = *(mat1->data);
    double * B = *(mat2->data);
    double * C = *(result->data);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m1r; i++) {
        for (int j = 0; j < m2c_e; j+=4*UNROLL) {
            __m256d c[UNROLL];

            for (int x = 0; x < UNROLL; x ++)
                c[x] = _mm256_setzero_pd();

            for (int k = 0; k < m1c; k++) {
                __m256d a = _mm256_broadcast_sd(A + i * m1c + k);
                for (int x = 0; x < UNROLL; x++) {
                    __m256d b = _mm256_loadu_pd(B + k * m2c + (j+4*x));
                    c[x] = _mm256_fmadd_pd(a, b, c[x]);
                }
            }
            
            for (int x = 0; x < UNROLL; x++)
                _mm256_storeu_pd(C + i * m2c + (j + x * 4), c[x]);
        }
    }    

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m1r; ++i) {
        for (int j = m2c_e; j < m2c; ++j) {
            double s = 0.0;
            for (int k = 0; k < m1c; ++k) {              
                s += A[i * m1c + k] * B[k * m2c + j];
            }
            C[i * m2c + j] = s;
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
    matrix ** tt = malloc(sizeof(matrix *));
    allocate_matrix(tt, mat->rows, mat->cols);
    matrix * t = *tt;

    matrix ** xx = malloc(sizeof(matrix *));
    allocate_matrix(xx, mat->rows, mat->cols);
    matrix * x = *xx;

    copy(mat, x);

    #pragma omp parallel for
    for(int i = 0; i < mat->rows; i++) {
        (*(result->data))[i + i * result->cols] = 1.0;
    }

    while (pow) {
        if (pow % 2 != 0) {
            mul_matrix_pow(t, result, x);
            copy(t, result);
            pow -= 1;
        }

        mul_matrix_pow(t, x, x);
        copy(t, x);
        pow /= 2;
    }

    deallocate_matrix(*tt);
    free(tt);

    deallocate_matrix(*xx);
    free(xx);
    
    return 0;
}


void copy(matrix * src, matrix* dest) {
    memcpy(*(dest->data), *(src->data), dest->rows * dest->cols * sizeof(double));
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    int rows = result->rows;
    int cols = result->cols;

    double *A = *(mat->data);
    double *B = *(result->data);

    int boundary = cols / 4 * 4;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<boundary; j+= 4) {
            int offset = i * cols + j;
            __m256d a = _mm256_loadu_pd(A + offset);
            _mm256_storeu_pd(B + offset, _mm256_xor_pd(a, _mm256_set1_pd(-0.0)));
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = boundary; j<cols; j++) {
            int offset = i * cols + j;
            B[offset] = -A[offset];
        }
    }
    return 0;

}


inline __m256d _m256_abs_pd(__m256d a) {
    const __m256d abs_mask = _mm256_set1_pd(-0.);
    return _mm256_andnot_pd(abs_mask, a); 
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    int rows = result->rows;
    int cols = result->cols;

    double *A = *(mat->data);
    double *B = *(result->data);

    int boundary = cols / 4 * 4;
    
    const __m256d b = _mm256_set1_pd(-1.0);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = 0; j<boundary; j+= 4) {
            int offset = i * cols + j;
            __m256d a = _mm256_loadu_pd(A + offset);
            _mm256_storeu_pd(B + offset, _mm256_max_pd(_mm256_mul_pd(b, a), a));
            // _mm256_storeu_pd(B + offset, _m256_abs_pd(a));
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i<rows; i++) {
        for (int j = boundary; j<cols; j++) {
            int offset = i * cols + j;
            B[offset] = fabs(A[offset]);
        }
    }
    return 0;
}

