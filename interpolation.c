//
//  interpolation.c
//  
//
//  Created by Laura Herrle on 12/7/15.
//
//

#include "interpolation.h"

#define NUM_ITER 100

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
double lagrange_interp_2d ( int mx, int my, double xd_1d[], double yd_1d[],
                            double zd[], double xi, double yi )

{
    int i;
    int j;
    int k;
    int l;
    double lx;
    double ly;
    double zi;
    
    
    l = 0;
    zi = 0.0;
    for ( j = 0; j < my + 1; j++ )
    {
        for ( i = 0; i < mx + 1; i++ )
        {
            lx = lagrange_basis_function_1d ( mx, xd_1d, i, xi );
            ly = lagrange_basis_function_1d ( my, yd_1d, j, yi );
            zi = zi + zd[l] * lx * ly;
            l = l + 1;
        }
    }
    
    return zi;
}

void matrix_gen (double x_0, double dx, int n_x, double v_0, double dv, int n_v, double dt, int n_l, double x[n_x], double v[n_v], double F[n_x+n_l*2][n_v+n_l*2]) {
    int i, j;
    for (i = 0; i < n_x; i++) {
        x[i+n_l] = x_0 + dx*i;
        for (j = 0; j < n_v; j++) {
            if (i == 0) {
                v[j+n_l] = v_0 + dv*j;
            }
            F[i+n_l][j+n_l] = cos(x[i+n_l]) * exp(-1*v[j+n_l]*v[j+n_l]);
            if(F[i+n_l][j+n_l] == 0) F[i+n_l][j+n_l] = 0;
        }
    }
    
    for (i = 0; i < n_l; i++) {
        x[i] = x_0 - dx*(n_l - i);
        x[n_l + n_x + i] = x_0 + dx*(n_x + i);
        for (j = 0; j < n_v + n_l*2; j++) {
            F[i][j] = F[n_x+i][j];
            F[i+n_x+n_l][j] = F[n_l+i][j];
        }
    }
    
    for (i = 0; i < n_l; i++) {
        v[i] = v_0 - dv*(n_l - i);
        v[n_l + n_v + i] = v_0 + dv*(n_v + i);
        for (j = 0; j < n_x + n_l*2; j++) {
            F[j][i] = 0;
            F[j][i+n_v+n_l] = 0;
        }
    }
}

void print_center_to_file(char* file, int n_x, int n_v, int n_l, double mat[n_x + n_l*2][n_v+n_l*2]) {
    FILE* fp;
    char* mode = "w";
    fp = fopen(file, mode);
    
    if(fp == NULL) {
        fprintf(stderr, "Can't open output file %s!\n",
                file);
        exit(1);
    }
    
    int i,j;
    
    for (i = 0; i < n_v; i ++) {
        for (j = 0; j < n_x; j ++) {
            if (j == n_x - 1) {
                fprintf(fp, "%.8f\n", mat[n_l + j][n_l+i]);
            } else {
            fprintf(fp, "%.8f, ", mat[n_l + j][n_l+i]);
            }
        }
    }
}

int main(int argc, char** argv) {
    double x_0 = 0;
    double dx = 0.1;
    int n_x = atoi(argv[1]);
    double v_0 = -.5;
    double dv = 0.01;
    int n_v = n_x;
    double dt = 0.0001;
    int n_l = atoi(argv[2]);
    double x[n_x+2*n_l];
    double v[n_v+2*n_l];
    double F_1[n_x+n_l*2][n_v+n_l*2];
    double F_2[n_x+n_l*2][n_v+n_l*2];
    double sub[n_l][n_l];
    double sub_x[n_l], sub_v[n_l];
    double new_point;
    
    matrix_gen(x_0, dx, n_x, v_0, dv, n_v, dt, n_l, x, v, F_1);
    
    print_center_to_file("F-init.csv", n_x, n_v, n_l, F_1);
    print_center_to_file("F-vs.csv", n_v+2*n_l, 1, 0, v);
    print_center_to_file("F-xs.csv", n_x+2*n_l, 1, 0, x);
    
    int i,j, x_c, v_c, ii, jj, iter;
    double x_tilde, v_tilde;
    
    for (iter = 0; iter < NUM_ITER; iter ++) {
        //on even iterations interpolate using f1 and put in f2
        #pragma omp parallel
        {
            if (iter%2==0) {
                //#pragma omp for
                for (i = n_l; i < n_x + n_l; i++) {
                    for (j = n_l; j < n_v + n_l; j++) {
                        x_tilde = x[i] + v[j]*dt;
                        v_tilde = v[j] + cos(x[i])*dt/2; //THIS IS WHERE E IS
                        x_c = floor((x_tilde-x_0)/dx) + n_l/2 + (n_l%2!=0);
                        v_c = floor((v_tilde-v_0)/dv) + n_l/2 + (n_l%2!=0);
                        
                        for (ii=0; ii<n_l; ii++) {
                            sub_x[ii] = x[x_c+ii];
                            for (jj = 0; jj < n_l; jj++) {
                                if(ii==0) {
                                    sub_v[jj] = v[v_c+jj];
                                }
                                sub[ii][jj] = F_1[x_c+ii][v_c+jj];
                            }
                        }
                        
                        new_point = lagrange_interp_2d(n_l-1, n_l-1, sub_x, sub_v, sub, x_tilde, v_tilde);
                        F_2[i][j] = new_point;
                    }
                }
            } else { //on odd do opposite
                //#pragma omp for
                for (i = n_l; i < n_x; i++) {
                    for (j = n_l; j < n_v; j++) {
                        x_tilde = x[i] + v[j]*dt;
                        v_tilde = v[j] + cos(x[i])*dt/2; //THIS IS WHERE E IS
                        x_c = floor((x_tilde-x_0)/dx) + n_l/2 + (n_l%2!=0);
                        v_c = floor((v_tilde-v_0)/dv) + n_l/2 + (n_l%2!=0);
                        
                        for (ii=0; ii<n_l; ii++) {
                            sub_x[ii] = x[x_c+ii];
                            for (jj = 0; jj < n_l; jj++) {
                                if(ii==0) {
                                    sub_v[jj] = v[v_c+jj];
                                }
                                sub[ii][jj] = F_2[x_c+ii][v_c+jj];
                            }
                        }
                        
                        new_point = lagrange_interp_2d(n_l-1, n_l-1, sub_x, sub_v, sub, x_tilde, v_tilde);
                        F_1[i][j] = new_point;
                    }
                }
            }
        }
    }
    
    if (NUM_ITER%2 == 0) {
        print_center_to_file("F-final.csv", n_x, n_v, n_l, F_2);
    } else {
        print_center_to_file("F-final.csv", n_x, n_v, n_l, F_1);
    }
}