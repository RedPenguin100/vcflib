/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/
#pragma once


double i4_binomial_pdf ( int n, double p, int k );
int i4_binomial_sample ( int n, double pp );
double i4vec_multinomial_pdf ( int n, double p[], int m, int x[] );
int *i4vec_multinomial_sample ( int n, double p[], int ncat );
double r8_beta_pdf ( double alpha, double beta, double rval );
double r8_beta_sample ( double aa, double bb );
double r8_chi_pdf ( double df, double rval );
double r8_chi_sample ( double df );
double r8_choose ( int n, int k );
double r8_epsilon ( void );
double r8_exponential_pdf ( double beta, double rval );
double r8_exponential_sample ( double lambda );
double r8_exponential_01_pdf ( double rval );
double r8_exponential_01_sample ( );
double r8_gamma_log ( double x );
double r8_gamma_pdf ( double beta, double alpha, double rval );
double r8_gamma_sample ( double a, double r );
double r8_gamma_01_pdf ( double alpha, double rval );
double r8_gamma_01_sample ( double a );
double r8_invchi_pdf ( double df, double rval );
double r8_invchi_sample ( double df );
double r8_invgam_pdf ( double beta, double alpha, double rval );
double r8_invgam_sample ( double beta, double alpha );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_normal_pdf ( double av, double sd, double rval );
double r8_normal_sample ( double av, double sd );
double r8_normal_01_pdf ( double rval );
double r8_normal_01_sample ( );
double r8_scinvchi_pdf ( double df, double s, double rval );
double r8_scinvchi_sample ( double df, double s );
double r8_uniform_pdf ( double lower, double upper, double rval );
double r8_uniform_sample ( double low, double high );
double r8_uniform_01_pdf ( double rval );
double r8_uniform_01_sample ( void );
double *r8mat_mtv_new ( int m, int n, double a[], double x[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
double r8mat_podet ( int n, double r[] );
double *r8mat_pofac ( int n, double a[] );
double *r8mat_poinv ( int n, double r[] );
double *r8mat_upsol ( int n, double r[], double b[] );
double *r8mat_utsol ( int n, double r[], double b[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_multinormal_pdf ( int n, double mu[], double r[], double c_det, 
  double x[] );
double *r8vec_multinormal_sample ( int n, double mu[], double r[] );

