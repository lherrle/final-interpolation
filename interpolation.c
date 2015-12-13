//
//  interpolation.c
//  
//
//  Created by Laura Herrle on 12/7/15.
//
//

#include "interpolation.h"

/*
 from lagrange_interp_2d.c by John Burkardt
 */
double lagrange_basis_function_1d ( int mx, double xd[], int i, double xi )
{
    int j;
    double yi;
    
    yi = 1.0;
    
    if ( xi != xd[i] )
    {
        for ( j = 0; j < mx + 1; j++ )
        {
            if ( j != i )
            {
                yi = yi * ( xi - xd[j] ) / ( xd[i] - xd[j] );
            }
        }
    }
    
    return yi;
}

/*
 from lagrange_interp_2d.c by John Burkardt
 */
double *lagrange_interp_2d ( int mx, int my, double xd_1d[], double yd_1d[],
                            double zd[], int ni, double xi[], double yi[] )

{
    int i;
    int j;
    int k;
    int l;
    double lx;
    double ly;
    double *zi;
    
    zi = ( double * ) malloc ( ni * sizeof ( double ) );
    
    for ( k = 0; k < ni; k++ )
    {
        l = 0;
        zi[k] = 0.0;
        for ( j = 0; j < my + 1; j++ )
        {
            for ( i = 0; i < mx + 1; i++ )
            {
                lx = lagrange_basis_function_1d ( mx, xd_1d, i, xi[k] );
                ly = lagrange_basis_function_1d ( my, yd_1d, j, yi[k] );
                zi[k] = zi[k] + zd[l] * lx * ly;
                l = l + 1;
            }
        }
    }
    return zi;
}

void matrix_gen (double x_0, double dx, int n_x, double v_0, double dv, int n_v, double dt, int n_l, double x[n_x], double v[n_v], double F[n_x][n_v]) {
    int i, j;
    for (i = 0; i < n_x; i++) {
        x[i] = x_0 + dx*i;
        for (j = 0; j < n_v; j++) {
            if (i == 0) {
                v[j] = v_0 + dv*j;
            }
            F[i][j] = cos(x[i]) * exp(-v[j]*v[j]);
        }
    }
}

int main(int argc, char** argv) {
    
}