#ifndef MATRIX_H
    #define MATRIX_H

    typedef struct matrix {
        int rows;
        int cols;
        double *data;
    } matrix;

    typedef struct regression_fit {
        int n;        
        int dof;

        double s2;
        double r2;
        double adjusted_r2;

        matrix *beta;
        matrix *beta_cov;
        matrix *beta_tstats;
    } regression_fit;


    // Indexing functions.
    double at(int, int, matrix*);
    double *ptr_at(int, int, matrix*);

    // Creation, allocation, and deletion functions
    matrix *new_matrix(int i, int j);
    matrix *copy_matrix(matrix*);
    matrix *make_matrix(int, int, ...);
    int init_matrix(double, matrix*);
    int free_data(matrix*);
    int free_matrix(matrix*);

    regression_fit *new_regression_fit(void);
    int free_regression_fit(regression_fit*);

    // Linear algebra functions.
    int transpose(matrix*, matrix*);
    int print_matrix(matrix*);
    int print_regression_fit(regression_fit*);
    int dot(matrix*, matrix*, matrix*);
    int matrix_add(matrix*, matrix*, matrix*);
    int matrix_subtract(matrix*, matrix*, matrix*);
    int matrix_scalar_multiply(matrix*, double, matrix*);
    int inverse(matrix*, matrix*);
    int lstsq(matrix*, matrix*, regression_fit*);

    double matrix_mean(matrix*, int);
    

#endif
