#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main(void){

    // Least squares test
    matrix *y = make_matrix(10, 1,
                             -0.05911,   
                              8.11098065,
                              2.24198061,
                             -8.0514162,
                             -8.07322677,
                             10.44922017,
                              8.82230702,
                             -3.18227219,
                              4.66138728, 
                            -21.56386711);

    matrix *x = make_matrix(10, 2,
                            1.0, -0.18380953,
                            1.0,  0.54366253,
                            1.0,  0.14948514,
                            1.0, -0.99259025,
                            1.0, -0.82037301,
                            1.0,  0.79391838,
                            1.0,  0.91688518,
                            1.0, -0.51775223,
                            1.0,  0.45836676,
                            1.0, -2.18044347);

    regression_fit *rfit = new_regression_fit();

    lstsq(y, x, rfit);
    print_regression_fit(rfit);

    free_matrix(y);
    free_matrix(x);
    free_regression_fit(rfit);
    return 0;
}
