#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "matrix.h"

/* TODO:
   1. Migrate shape variables from int to size_t or unsigned.
   2. Make a helper function for null checking.
   3. Add checks for shape conformity.
   4. Add checks in Gaussian elimination for singularity.
   5. Better way to share and/or free matrices with regression_fit? 
 */

/************************************
 * Private implementation functions *
 ************************************/

static int 
index(int i, int j, int row_length){
    return (i * row_length) + j;
}

static int 
m_init_matrix(double init_val, matrix *m){
    for(int i=0; i < m->rows * m->cols; i++){
        m->data[i] = init_val;
    }
    return 1;
}

static int 
m_free_data(matrix *m){
    free(m->data);
    return 1;
}

static int
m_free_matrix(matrix *m){
    free(m);
    return 1;
}

static int 
m_print_matrix(matrix *m){
    printf("[%d x %d]\n[", m->rows, m->cols);
    for(int i=0; i < m->rows; i++){
        
        // Prepend space for left alignment.
        if(i > 0){
            printf(" ");
        }
        
        // Print each entry for a given row.
        for(int j=0; j < m->cols; j++){
            printf("%.5f", at(i, j, m));
            
            // If at the end of data, print a closing bracket.
            if(i == m->rows-1 && j == m->cols-1){
                printf("]");   
            } 
            // Otherwise print a comma delimiter.
            else {
                printf(", ");
            }
        }

        // End every row with a newline.
        printf("\n");
    }

    return 1;
}

static int 
m_print_regression_fit(regression_fit *rfit){


    printf("-------------------------"
           "Summary of Regression Analysis"
           "-------------------------\n\n");

    printf("Number of Observations:         %d\n", rfit->n);
    printf("Number of Degrees of Freedom:   %d\n\n", rfit->dof);
    printf("R-squared:         %.4f\n", rfit->r2);
    printf("Adj R-squared:     %.4f\n\n", rfit->adjusted_r2);

    printf("-----------------------"
           "Summary of Estimated Coefficients"
           "------------------------\n");

    printf("     Variable       Coef      t-stat     \n");

    printf("---------------------------"
           "---------------------------"
           "--------------------------\n");


    for(int i=0; i < rfit->beta->rows; i++){
        printf("     Column_%d     %.4f      %.4f     \n", 
               i,
               at(i, 0, rfit->beta),
               at(i, 0, rfit->beta_tstats));            
    }

    printf("---------------------------------"
           "End of Summary"
           "---------------------------------\n\n");

    return 1;
}

static int 
m_transpose(matrix *m, matrix *m_t){
    m_t->rows = m->cols;
    m_t->cols = m->rows;

    int position = 0;
    for(int j=0; j < m->cols; j++){
        for(int i=0; i < m->rows; i++){
            m_t->data[position] = at(i, j, m);
            position++;
        }
    }

    return 1;
}

static int 
m_dot(matrix *a, matrix *b, matrix *c){
    double val = 0.0;
    
    for(int i=0; i < a->rows; i++){
        for(int k=0; k < b->cols; k++){
            val = 0.0;

            for(int j=0; j < b->rows; j++){
                val += at(i, j, a) * at(j, k, b);
            }

            *ptr_at(i, k, c) = val;
        }
    }

    return 1;
}

static int 
m_matrix_add(matrix *a, matrix *b, matrix *c){
    for(int i=0; i < a->rows; i++){
        for(int j=0; j < a->cols; j++){
            *ptr_at(i, j, c) = at(i, j, a) + at(i, j, b);
        }
    }

    return 1;
}

static int 
m_matrix_subtract(matrix *a, matrix *b, matrix *c){
    for(int i=0; i < a->rows; i++){
        for(int j=0; j < a->cols; j++){
            *ptr_at(i, j, c) = at(i, j, a) - at(i, j, b);
        }
    }

    return 1;
}

static int 
m_matrix_scalar_multiply(matrix *a, double b, matrix *c){
    for(int i=0; i < a->rows; i++){
        for(int j=0; j < a->cols; j++){
            *ptr_at(i, j, c) = b * at(i, j, a);
        }
    }

    return 1;
}

static double 
m_matrix_mean(matrix *a, int col){
    double result = 0.0;

    /* Tabulate the sum along column 'col'.*/
    for(int i=0; i < a->rows; i++){
        result += at(i, col, a);
    }

    return result / ( (double) a->rows);
}

static int 
m_gauss_solve(matrix *a_in, matrix *b, matrix *x){
    int i; /* Used for loops over rows and columns. */
    int j; /* Used for loops over columns. */
    int k; /* Used for loops over rows. */
    int m; /* Used for identifying the pivot row at each step. */

    double row_multiplier; /* Multiple of row to be subtracted at each step. */
    double temp;           /* Placeholder for swap/pivot operations. */ 
    double pivot_value;    /* Container for max pivot value. */

    int n = a_in->rows;

    /* Allocate and initialize a copy, so that the process of
     * solving Ax=b doesn't mutate A's data. Be sure to free. */
    matrix *a = copy_matrix(a_in);

    /* Forward steps to pivot rows and perform subtraction of row
     * multiples. 
     *
     * k tracks columns in forward order.
     * m is used to track the row with maximal pivot entry in kth column.
     * i is used to track rows that follow the kth row.
     * j is used to track columns that follow the kth column. */
    for (k=0; k < n-1; k++) {
	
        pivot_value = (double) fabs(at(k, k, a));
        m = k; /* m is pivot row, assumed equal to current row. */

        /* Find the remaining row with largest pivot */
        for (i=k+1; i < n; i++){   
            temp = (double) fabs(at(i, k, a));            
            if(temp > pivot_value){
                pivot_value = temp; 
                m = i; /* m is set to row with largest pivot in kth column. */
            }
        }

        // If location, m, of largest pivot element is not row k, then
        // swap row m and row k in a and b.
        if(m != k) {  

            /* Swap entries of b vector. */
            temp = at(k, 0, b);
            *ptr_at(k, 0, b)  = at(m, 0, b);
            *ptr_at(m, 0, b)  = temp;

            /* Here j is used for swapping the columns between row k and row m.
             * This operation can effectively ignore swapping entries in the
             * lower triangle, which is why j goes from k to n-1. The nature
             * of the Gaussian elimination algorithm is such that by the end,
             * all lower triangle elements must be zero, and are ignored. */

            for(j=k; j < n; j++) {
                temp = at(k, j, a);
                *ptr_at(k, j, a) = at(m, j, a);
                *ptr_at(m, j, a) = temp;
            }
        }

        /* For each remaining row after the kth row, we are going to subtract a
         * multiple of the pivot row from that row. The only columns affected 
         * are columns (k+1) through the final column. */

        for (i=k+1; i < n; i++) {
            row_multiplier = at(i, k, a) / at(k, k, a);
            
            /* Notice how the loop below only affects columns to the right of
             * the current column (k). This is because the formula would 
             * implicitly zero-out the entry at (j,k) -- the multiplier is
             * a[i, k] / a[k, k] -- and when multiplied by a[k, k] it leaves
             * just a[i, k]. So when j=k, the subtraction would result in 
             * a[i, k] - a[i, k] = 0. 
             *
             * Rather than doing this redundant computation. The values in the
             * lower triangle of a are simply ignored and never used during the
             * back-substitution phase, treating them as 0 without computing
             * with them. */

            for (j=k+1; j < n; j++) {
                *ptr_at(i, j, a) = at(i, j, a) - row_multiplier * at(k, j, a);
            }

            /* Adjust the RHS b vector by the same operation.*/
            *ptr_at(i, 0, b) = at(i, 0, b) - row_multiplier * at(k, 0, b);
        }       
    }

    /* Back substitution to get final solution.
     * j tracks the number of rows.
     * k is used to iterate the rows in reverse.
     * i is used to iterate upper-triangle columns of a given row. */

    for (j=0; j < n; j++) {

        /* k tracks from the final row back to the first row. */
        k = n - j - 1;

        /* Begin by setting the value of x[k] to b[k]. */
        *ptr_at(k, 0, x) = at(k, 0, b);

        /* By the time we get to column i, the value of x is already solved for
         * row i, so the value x[i] is available to use. The scalar multiple of
         * x[i] in the equation we want to solve for x[k] will be a[k, i], so we
         * subtract a[k,i] * x[i] from the running solution for x[k] */

        for(i=k+1; i < n; i++) {
            *ptr_at(k, 0, x) = at(k, 0, x) - at(k, i, a) * at(i, 0, x);
        }

        /* Finally, we divide by the pivot entry in row k, since everything
         * else has correctly been accounted for in x[k]. */

        *ptr_at(k, 0, x) = at(k, 0, x) / at(k, k, a);

        /* Note that the above loop would be skipped when k=n-1, for example 
         * and x[n-1] would just be set to the value b[n-1]/a[n-1, n-1] at that
         * point. Then the loop would only have one iteration for k=n-2, and
         * x[n-2] would be set to: 
         *     (b[n-2] - a[n-2, n-1] * x[n-1]) / a[n-2, n-2] 
         * ... and so forth back up to k=0. */
    }

    free_matrix(a);
    return 1;
}

static int 
m_inverse(matrix *m, matrix *m_inv){
    
    matrix *x = new_matrix(m->rows, 1);
    matrix *b = new_matrix(m->rows, 1);

    // Solve Ax=b once for each column.
    for(int c=0; c < m->cols; c++){
        
        // Initialize to c-th coordinate vector,
        // in case pivoting changed it.
        for(int i=0; i < m->rows; i++){
            if(i == c){
                *ptr_at(i, 0, b) = 1.0;
            } else {
                *ptr_at(i, 0, b) = 0.0;
            }
        }

        m_gauss_solve(m, b, x);

        for(int i=0; i < m->rows; i++){
            *ptr_at(i, c, m_inv) = at(i, 0, x);
        }
    }

    free_matrix(x);
    free_matrix(b);

    return 1;
}

static int 
m_lstsq(matrix *y, matrix *x, regression_fit *rfit){

    matrix *xT      = new_matrix(x->cols, x->rows);
    matrix *xTy     = new_matrix(x->cols, y->cols);
    matrix *xTx     = new_matrix(x->cols, x->cols);
    matrix *xTx_inv = new_matrix(x->cols, x->cols);

    matrix *eps     = new_matrix(y->rows, y->cols);
    matrix *epsT    = new_matrix(y->cols, y->rows);

    matrix *y_hat   = new_matrix(y->rows, y->cols);
    matrix *epsTeps = new_matrix(y->cols, y->cols);

    int n = x->rows;
    int degrees_of_freedom = x->cols;
    double y_bar = matrix_mean(y, 0);

    double ssr = 0.0;
    double tss = 0.0;
    double n_m_1;
    double n_m_p_1;
    double std_err;

    rfit->n   = n;
    rfit->dof = degrees_of_freedom;

    // For regression fit.
    rfit->beta = new_matrix(x->cols, y->cols);

    transpose(x, xT);
    dot(xT, x, xTx);
    inverse(xTx, xTx_inv);
    dot(xT, y, xTy);
    dot(xTx_inv, xTy, rfit->beta);

    // For summary statistics.
    dot(x, rfit->beta, y_hat);
    matrix_subtract(y, y_hat, eps);
    transpose(eps, epsT);
    dot(epsT, eps, epsTeps);

    rfit->s2 = at(0, 0, epsTeps) / (n - degrees_of_freedom);

    // Sum of squares for r2 calculation.
    for(int i=0; i < y->rows; i++){
        ssr += pow(at(i, 0, y_hat) - y_bar, 2.0);
        tss += pow(at(i, 0, y) - y_bar, 2.0);
    }
    
    rfit->r2 = (ssr / tss);

    // Adjusted R2
    n_m_1   = rfit->n - 1;
    n_m_p_1 = rfit->n - (rfit->dof - 1) - 1; 
    rfit->adjusted_r2 = 1.0 - (1.0 - rfit->r2) * (n_m_1 / n_m_p_1); 
    
    // Beta covariance and t-stats
    rfit->beta_cov = new_matrix(xTx_inv->rows, xTx_inv->cols);
    matrix_scalar_multiply(xTx_inv, rfit->s2, rfit->beta_cov);

    rfit->beta_tstats = new_matrix(rfit->beta->rows, rfit->beta->cols);
    for(int i=0; i < rfit->beta->rows; i++){
        std_err = sqrt(at(i, i, rfit->beta_cov));
        *ptr_at(i, 0, rfit->beta_tstats) = at(i, 0, rfit->beta) / std_err;
    }

    free_matrix(xT);
    free_matrix(xTy);
    free_matrix(xTx);
    free_matrix(xTx_inv);
    free_matrix(eps);
    free_matrix(epsT);
    free_matrix(epsTeps);
    free_matrix(y_hat);
    
    return 1;
}



/************************
 * Public API Functions *
 ************************/

double 
at(int i, int j, matrix *m){
    if(m == NULL){
        return 0.0;
    } else if(m->data == NULL){
        return 0.0;
    } else if(i >= m->rows || j >= m->cols){
        return 0.0;
    }

    return m->data[index(i, j, m->cols)];
}

double *
ptr_at(int i, int j, matrix *m){
    if(m == NULL){
        return NULL;
    } else if(m->data == NULL){
        return NULL;
    } else if(i >= m->rows || j >= m->cols){
        return NULL;
    }

    return &(m->data[index(i, j, m->cols)]);
}

int 
init_matrix(double init_val, matrix *m){
    if(m == NULL){
        return -1;
    } else if(m->data == NULL){
        return -2;
    }
 
    return m_init_matrix(init_val, m);
}

matrix *
new_matrix(int i, int j){
    matrix *m = (matrix *) malloc(sizeof(matrix));
    if(m == NULL){
        return NULL;
    } 
    
    m->rows = i;
    m->cols = j;
    m->data = (double *) malloc(sizeof(double) * m->rows * m->cols);
    
    if(m->data == NULL){
        free(m);
        return NULL;
    }

    m_init_matrix(0.0, m);
    return m;
}

matrix *
copy_matrix(matrix *m){
    if(m == NULL){
        return NULL;
    }

    if(m->data == NULL){
        return NULL;
    }

    matrix *output = new_matrix(m->rows, m->cols);
    for(int i=0; i < m->rows * m->cols; i++){
        output->data[i] = m->data[i];
    }

    return output;
}

matrix *
make_matrix(int i, int j, ...){
    va_list arg_list;
    va_start(arg_list, j);

    matrix *m = (matrix *) malloc(sizeof(matrix));
    if(m == NULL){
        return NULL;
    }
    
    m->data = (double *) malloc(sizeof(double) * i * j);
    
    if(m->data == NULL){
        free(m);
        return NULL;
    }

    m->rows = i;
    m->cols = j;

    for(int pos=0; pos < i*j; pos++){
        m->data[pos] = va_arg(arg_list, double);    
    }

    va_end(arg_list);
    return m;
}

int 
free_data(matrix *m){
    if(m == NULL){
        return -1;
    } else if(m->data == NULL){
        return -2;
    }
 
    return m_free_data(m);
}

int 
free_matrix(matrix *m){
    int data_status = free_data(m);
    if(data_status == 0){
        return m_free_matrix(m);
    } 

    return data_status;
}

regression_fit *
new_regression_fit(void){
    regression_fit *rfit = (regression_fit *) malloc(sizeof(regression_fit));
    if(rfit == NULL){
        return NULL;
    } 
    
    return rfit;
}

int 
free_regression_fit(regression_fit *rfit){
    if(   rfit->beta == NULL
       || rfit->beta_cov == NULL 
       || rfit->beta_tstats == NULL){
        return -1;
    }

    // TODO -- further check return statuses of free_matrix.
    free_matrix(rfit->beta);
    free_matrix(rfit->beta_cov);
    free_matrix(rfit->beta_tstats);
    free(rfit);

    return 1;
}

int 
print_matrix(matrix *m){
    if(m == NULL){
        return -1;
    } else if(m->data == NULL){
        return -2;
    }
    return m_print_matrix(m);
}

int 
print_regression_fit(regression_fit *rfit){
    if(rfit == NULL){
        return -1;
    }

    return m_print_regression_fit(rfit);
}

int 
transpose(matrix *m, matrix *m_t){
    if(m == NULL || m_t == NULL){
        return -1;
    } else if(m->data == NULL || m_t->data == NULL){
        return -2;
    } 
    
    return m_transpose(m, m_t);
}

int 
dot(matrix *a, matrix *b, matrix *c){
    if(a == NULL || b == NULL || c == NULL){
        return -1;
    } else if(a->data == NULL || b->data == NULL || c->data == NULL){
        return -2;
    } else if(a->cols != b->rows){
        return -3;
    } else if(a->rows != c->rows || b->cols != c->cols){
        return -4;
    }

    return m_dot(a, b, c);
}

int 
matrix_add(matrix *a, matrix *b, matrix *c){
    if(a == NULL || b == NULL || c == NULL){
        return -1;
    } else if(a->data == NULL || b->data == NULL || c->data == NULL){
        return -2;
    } else if(a->rows != b->rows || b->rows != c->rows){
        return -3;
    } else if(a->cols != b->cols || b->cols != c->cols){
        return -4;
    }

    return m_matrix_add(a, b, c);
}

int 
matrix_subtract(matrix *a, matrix *b, matrix *c){
    if(a == NULL || b == NULL || c == NULL){
        return -1;
    } else if(a->data == NULL || b->data == NULL || c->data == NULL){
        return -2;
    } else if(a->rows != b->rows || b->rows != c->rows){
        return -3;
    } else if(a->cols != b->cols || b->cols != c->cols){
        return -4;
    }

    return m_matrix_subtract(a, b, c);
}

int 
matrix_scalar_multiply(matrix *a, double b, matrix *c){
    if(a == NULL || c == NULL){
        return -1;
    } else if(a->data == NULL || c->data == NULL){
        return -2;
    } else if(a->rows != c->rows || a->cols != c->cols){
        return -3;
    }

    return m_matrix_scalar_multiply(a, b, c);
}

int 
inverse(matrix *a, matrix *a_inv){
    if(a == NULL || a_inv == NULL){
        return -1;
    } else if(a->data == NULL || a_inv->data == NULL){
        return -2;
    } else if(a->rows != a_inv->rows){
        return -3;
    } else if(a->cols != a_inv->cols){
        return -4;
    }

    return m_inverse(a, a_inv);
}

int 
lstsq(matrix *y, matrix *x, regression_fit *rfit){
    if(y == NULL || x == NULL || rfit == NULL){
        return -1;
    }

    if(y->data == NULL || x->data == NULL){
        return -2;
    }

    // TODO
    // add shape conformity and missing data checks.

    return m_lstsq(y, x, rfit);
}

double 
matrix_mean(matrix *a, int col){
    if(a == NULL || a->data == NULL){
        return 0.0;
    } else if(a->cols >= col){
        return 0.0;
    }

    return m_matrix_mean(a, col);
}


