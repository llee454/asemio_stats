#include <stdio.h>
#include <string.h> // for memcpy

// Prevent OCaml from exporting short macro names.
#define CAML_NAME_SPACE 1

#include <caml/callback.h>
#include <caml/mlvalues.h>
#include <caml/memory.h> // CAMLreturn
#include <caml/alloc.h> // caml_copy

#include <gsl_errno.h>
#include <gsl_math.h> // log1p
#include <gsl_statistics_double.h>
#include <gsl_sf_erf.h>
#include <gsl_fit.h> // gsl_fit_linear
#include <gsl_integration.h>
#include <gsl_multifit_nlinear.h>
#include <gsl_cdf.h>
#include <gsl_sf_gamma.h> // gsl_sf_fact
#include <gsl_randist.h> // gsl_ran_binomial_pdf

#include <ocaml_siman.h>

CAMLprim value ocaml_gsl_pow_int (value x, value n) {
  CAMLparam2 (x, n);
  CAMLreturn (caml_copy_double (gsl_pow_int (Double_val (x), Int_val (n))));
}

CAMLprim value ocaml_mean (value xs) {
  CAMLparam1 (xs);
  size_t n = Wosize_val (xs);
  double* x = malloc (n * sizeof (double));
  for (size_t i = 0; i < n; i ++) {
    x [i] = Double_field (xs, i);
  }
  CAMLreturn (caml_copy_double (gsl_stats_mean (x, 1, n)));
}

CAMLprim value ocaml_gsl_sf_erf_Z (value x) {
  CAMLparam1 (x);
  CAMLreturn (caml_copy_double (gsl_sf_erf_Z (Double_val (x))));
}

CAMLprim value ocaml_gsl_sf_erf_Q (value x) {
  CAMLparam1 (x);
  CAMLreturn (caml_copy_double (gsl_sf_erf_Q (Double_val (x))));
}

CAMLprim value ocaml_gsl_sf_fact (value x) {
  CAMLparam1 (x);
  CAMLreturn (caml_copy_double (gsl_sf_fact ((unsigned int) (Int_val (x)))));
}

CAMLprim value ocaml_gsl_cdf_gaussian_P (value x, value std) {
  CAMLparam2 (x, std);
  CAMLreturn (caml_copy_double (gsl_cdf_gaussian_P (Double_val (x), Double_val (std))));
}

CAMLprim value ocaml_gsl_cdf_gaussian_Q (value x, value std) {
  CAMLparam2 (x, std);
  CAMLreturn (caml_copy_double (gsl_cdf_gaussian_Q (Double_val (x), Double_val (std))));
}

CAMLprim value ocaml_gsl_cdf_chisq_P (value x, value nu) {
  CAMLparam2 (x, nu);
  CAMLreturn (caml_copy_double (gsl_cdf_chisq_P (Double_val (x), Double_val (nu))));
}

CAMLprim value ocaml_gsl_cdf_gaussian_Pinv (value x, value ndf) {
  CAMLparam2 (x, ndf);
  CAMLreturn (caml_copy_double (gsl_cdf_gaussian_Pinv (Double_val (x), Double_val (ndf))));
}

CAMLprim value ocaml_gsl_ran_binomial_pdf (value k, value p, value n) {
  CAMLparam3 (k, p, n);
  CAMLlocal1 (result);
  CAMLreturn (caml_copy_double (gsl_ran_binomial_pdf (
    (unsigned int) (Int_val (k)),
    Double_val (p),
    (unsigned int) (Int_val (n)))));
}

CAMLprim value ocaml_gsl_fit_linear (value xs, value ys) {
  CAMLparam2 (xs, ys);
  CAMLlocal1 (result);
  // BUG: wrong allocation size (Wosize * sizeof double)
  double* x = malloc (Wosize_val (xs));
  for (size_t i = 0; i < Wosize_val (xs); i ++) {
    x [i] = Double_field (xs, i);
  }
  double* y = malloc (Wosize_val (ys));
  for (size_t i = 0; i < Wosize_val (ys); i ++) {
    y [i] = Double_field (ys, i);
  }
  const size_t xstride = 1;
  const size_t ystride = 1;
  const size_t n = 3;
  double c0;
  double c1;
  double cov00;
  double cov01;
  double cov11;
  double sumsq;
  int status = gsl_fit_linear (x, xstride, y, ystride, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  free (x);
  free (y);

  // Note: OCaml represents records whose fields are all floats as double array blocks.
  result = caml_alloc (2, Double_array_tag);
  Store_double_field (result, 0, c0);
  Store_double_field (result, 1, c1);
  CAMLreturn (result);
}

struct callback_params { value h; };

double callback (double x, void* params) {
  CAMLparam0 ();
  CAMLlocal2 (y, result);
  struct callback_params* p = (struct callback_params*) params;
  y = caml_copy_double (x);
  result = caml_callback (p->h, y);
  return Double_val (result);
}

CAMLprim value ocaml_integrate (value f, value lower, value upper) {
  CAMLparam3 (f, lower, upper);
  CAMLlocal1 (result);
  double out;
  double err;
  size_t neval;
  struct callback_params params = {
    .h = f
  };
  gsl_function F = {
    .function = &callback,
    .params   = &params
  };
  int status = gsl_integration_qng (&F, Double_val (lower), Double_val (upper), 0.0001, 0, &out, &err, &neval);
  result = caml_alloc (3, 0);
  Store_field (result, 0, caml_copy_double (out));
  Store_field (result, 1, caml_copy_double (err));
  Store_field (result, 2, Val_long (neval));
  CAMLreturn (result);
}

CAMLprim value ocaml_integration_qag (value f, value lower, value upper) {
  CAMLparam3 (f, lower, upper);
  CAMLlocal1 (result);
  double out;
  double err;
  size_t limit = 10;
  struct callback_params params = {
    .h = f
  };
  gsl_function F = {
    .function = &callback,
    .params   = &params
  };
  gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);
  int status =  gsl_integration_qag (&F, Double_val (lower), Double_val (upper), 0.0001, 0, limit, GSL_INTEG_GAUSS61, w, &out, &err);
  gsl_integration_workspace_free (w);
  result = caml_alloc (3, 0);
  Store_field (result, 0, caml_copy_double (out));
  Store_field (result, 1, caml_copy_double (err));
  Store_field (result, 2, Val_long (0));
  CAMLreturn (result);
}

CAMLprim value ocaml_integration_qagi (value f) {
  CAMLparam1 (f);
  CAMLlocal1 (result);
  double out;
  double err;
  size_t limit = 10;
  struct callback_params params = {
    .h = f
  };
  gsl_function F = {
    .function = &callback,
    .params   = &params
  };
  gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);
  int status =  gsl_integration_qagi (&F, 0.0001, 0, limit, w, &out, &err);
  result = caml_alloc (3, 0);
  Store_field (result, 0, caml_copy_double (out));
  Store_field (result, 1, caml_copy_double (err));
  Store_field (result, 2, Val_long (0));
  gsl_integration_workspace_free (w);
  CAMLreturn (result);
}

// an auxiliary data structure used to store the parameters and data points fitted by ocaml_get_fit_nlinear
struct ocaml_get_fit_nlinear_f_data {
  size_t  n_xs;
  size_t  n_ks;
  value   f;
  double* xs;
  double* ys;
};

// computes the value of the fitted function for the given function parameters and data values.
int ocaml_get_fit_nlinear_f (
  const gsl_vector* f_ks, // the values of constants within the fitted function
  void* data,             // the fitted function parameters and data values
  gsl_vector* f_res       // the fitted function output values (may be 2 or more if f is multidimensional)
) {
  CAMLparam0 ();
  CAMLlocal3 (ocaml_f, ocaml_f_ks, ocaml_f_arg);

  size_t n_xs = ((struct ocaml_get_fit_nlinear_f_data*) data)->n_xs;   // the number of data points
  size_t n_ks = ((struct ocaml_get_fit_nlinear_f_data*) data)->n_ks;   // the number of fitted function parameters
  ocaml_f     = ((struct ocaml_get_fit_nlinear_f_data*) data)->f; // the fitted function
  double* xs  = ((struct ocaml_get_fit_nlinear_f_data*) data)->xs;
  double* ys  = ((struct ocaml_get_fit_nlinear_f_data*) data)->ys;

  // the fitted function parameters
  ocaml_f_ks = caml_alloc (n_ks, Double_array_tag);
  for (size_t i = 0; i < 2; i ++) {
    Store_double_field (ocaml_f_ks, i, gsl_vector_get (f_ks, i));
  }

  // call the fitted function on the given data values with the current parameter values.
  for (size_t i = 0; i < n_xs; i ++) {
    // construct the argument passed to the fitted function
    ocaml_f_arg = caml_alloc (2, 0);
    Store_field (ocaml_f_arg, 0, ocaml_f_ks);
    Store_field (ocaml_f_arg, 1, caml_copy_double (xs [i]));

    // call the fitted function and store the result
    gsl_vector_set (f_res, i, ys [i] - Double_val (caml_callback (ocaml_f, ocaml_f_arg)));
  }

  return GSL_SUCCESS;
}

/*
  Accepts a function that has parameters ks, and a dataset ys over over
  values xs and estimates the parameter values ks' that makes f
  approximate ys over xs.
*/
CAMLprim value ocaml_gsl_fit_nlinear (
  value ocaml_f,
  value ocaml_f_ks_init,
  value ocaml_xs,
  value ocaml_ys
) {
  CAMLparam4 (ocaml_f, ocaml_f_ks_init, ocaml_xs, ocaml_ys);
  CAMLlocal2 (ocaml_f_ks_final, ocaml_f_arg);

  const size_t max_iter = 100; // maximum number of iterations
  const double x_tol = 0.00001; // x error tolerance
  const double g_tol = 0.00001; // error tolerance
  const double f_tol = 0.00001; // f error tolerance
  int info; // the reason for convergence if curve fitting succeeds

  size_t n_xs = Wosize_val (ocaml_xs);
  size_t n_ys = Wosize_val (ocaml_ys);
  size_t n = n_xs < n_ys ? n_xs : n_ys;
  size_t p = Wosize_val (ocaml_f_ks_init);

  double* xs = malloc (n * sizeof (double));
  double* ys = malloc (n * sizeof (double));
  for (size_t i = 0; i < n; i ++) {
    xs [i] = Double_field (ocaml_xs, i);
    ys [i] = Double_field (ocaml_ys, i);
  }

  struct ocaml_get_fit_nlinear_f_data f_params = {
    .n_xs = n,
    .n_ks = p,
    .f    = ocaml_f,
    .xs   = xs,
    .ys   = ys
  };

  gsl_vector* ks_init = gsl_vector_alloc (p);
  for (size_t i = 0; i < p; i ++) {
    gsl_vector_set (ks_init, i, Double_field (ocaml_f_ks_init, i));
  }

  // configure the solver.
  const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;

  gsl_multifit_nlinear_fdf fdf;
  fdf.f      = ocaml_get_fit_nlinear_f;
  fdf.df     = NULL; // use finite-difference Jacobin (slope approximation)
  fdf.fvv    = NULL; // disable geodesic acceleration
  fdf.n      = n;
  fdf.p      = p;
  fdf.params = &f_params;

  // allocate workspace with initial solver parameters.
  gsl_multifit_nlinear_workspace* w;

  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  // initialize the solver.
  gsl_multifit_nlinear_init (ks_init, &fdf, w);

  int status = gsl_multifit_nlinear_driver (max_iter,
    x_tol, g_tol, f_tol, 
    NULL, // no callback function
    NULL, // no callback function params
    &info, w
  );

  // get the estimate fitted function parameters.
  ocaml_f_ks_final = caml_alloc (p, Double_array_tag);
  for (size_t i = 0; i < p; i ++) {
    Store_double_field (ocaml_f_ks_final, i,
      gsl_vector_get (w->x, i));
  }
 
  // clean up allocated memory resources.
  gsl_multifit_nlinear_free (w);
  gsl_vector_free (ks_init);
  free (xs);
  free (ys);

  CAMLreturn (ocaml_f_ks_final);
}

/*
struct ocaml_siman_callbacks {
  value copy_fn;
  value energy_fn;
  value step_fn;
  value distance_fn;
};

struct ocaml_siman_state {
  struct ocaml_siman_callbacks* callbacks;
  value state;
};

double ocaml_siman_energy (void* xp) {
  CAMLparam0 ();
  printf ("[ocaml_siman_energy]\n");
  struct ocaml_siman_state* s = (struct ocaml_siman_state*) xp;
  double result = Double_val (caml_callback (s->callbacks->energy_fn, s->state));
  printf ("[energy] x: %0.4f energy: %0.4f\n", Double_val (Field (s->state, 0)), result);
  fflush (stdout);
  CAMLreturnT (double, result);
}

void ocaml_siman_step (const gsl_rng* rng, void* xp, double step_size) {
  CAMLparam0 ();
  printf ("[ocaml_siman_step]\n");
  struct ocaml_siman_state* s = (struct ocaml_siman_state*) xp;
  double step   = (2 * gsl_rng_uniform (rng) - 1) * step_size;
  double orig_state = Double_val (Field (s->state, 0));
  caml_callback2 (s->callbacks->step_fn, s->state, caml_copy_double (step));
  printf (
    "[ocaml_siman_step] step: %0.4f orig: %0.4f next: %0.4f\n",
    step, orig_state, Double_val (Field (s->state, 0)));
  fflush (stdout);
  CAMLreturn0;
}

double ocaml_siman_distance (void* xp, void* yp) {
  CAMLparam0 ();
  printf ("[ocaml_siman_distance]\n");
  struct ocaml_siman_state* x = (struct ocaml_siman_state*) xp;
  struct ocaml_siman_state* y = (struct ocaml_siman_state*) yp;
  if (x->callbacks != y->callbacks) {
    GSL_ERROR("Internal Error in ocaml_siman", GSL_EINVAL);
  }
  double distance = Double_val (caml_callback2 (x->callbacks->distance_fn, x->state, y->state));
  printf ("[distance] result: %f\n", distance);
  fflush (stdout);
  CAMLreturnT (double, distance);
}

void ocaml_siman_copy (void* source, void* dest) {
  CAMLparam0 ();
  CAMLlocal1 (state);
  printf ("[ocaml_siman_copy]\n");
  fflush (stdout);
  struct ocaml_siman_state* s = (struct ocaml_siman_state*) source;
  struct ocaml_siman_state* d = (struct ocaml_siman_state*) dest;
  state = caml_callback (s->callbacks->copy_fn, s->state);
//   caml_register_global_root (state);
  d->callbacks   = (*s).callbacks;
  d->state       = state;
  printf ("[ocaml_siman_copy] done\n");
  fflush (stdout);
  CAMLreturn0 ();
}

void* ocaml_siman_construct (void* xp) {
  CAMLparam0 ();
  printf ("[ocaml_siman_construct]\n");
  fflush (stdout);
  struct ocaml_siman_state* dest = malloc (sizeof (struct ocaml_siman_state));
  ocaml_siman_copy (xp, dest);
  printf ("[ocaml_siman_construct] done\n");
  fflush (stdout);
  CAMLreturnT (struct ocaml_siman_state*, dest);
}

void ocaml_siman_destroy (void* xp) {
  CAMLparam0 ();
  CAMLlocal1 (state);
  printf ("[ocaml_siman_destroy]\n");
  fflush (stdout);
  struct ocaml_siman_state* x = (struct ocaml_siman_state*) xp;
  state = x->state;
//   caml_remove_global_root ((value*) &state);
//   free (xp);
  printf ("[ocaml_siman_destroy] done\n");
  fflush (stdout);
  CAMLreturn0;
}

void print_cfg (void* xp) {
  CAMLparam0 ();
  struct ocaml_siman_state* s = (struct ocaml_siman_state*) xp;
  printf (" print cfg x: %0.4f", Double_val (Field (s->state, 0)));
  fflush (stdout);
  CAMLreturn0;
} 

CAMLprim value ocaml_siman_solve (value copy_fn, value energy_fn, value step_fn, value distance_fn, value initial) {
  CAMLparam5 (copy_fn, energy_fn, step_fn, distance_fn, initial);

  printf ("[ocaml_siman_solve] initial: %0.4f\n", Double_val (Field (initial, 0)));
  fflush (stdout);

  // create the initial configuration state.
  struct ocaml_siman_callbacks callbacks = {
    .copy_fn     = copy_fn,
    .energy_fn   = energy_fn,
    .step_fn     = step_fn,
    .distance_fn = distance_fn
  };
  struct ocaml_siman_state xp0 = {
    .callbacks = &callbacks,
    .state     = initial
  };

  // printf ("[ocaml_siman_solve] callback s: %lu\n", (unsigned long) xp0.callbacks);

  // initialize a random number generator.
  gsl_rng_env_setup ();
  gsl_rng* rng = gsl_rng_alloc (gsl_rng_default);

  // paramters:
  gsl_siman_params_t params = {
    .n_tries = 2, // 200
    .iters_fixed_T = 8, // 1000
    .step_size = 1,
    .k = 1.0, // Boltzman constant
    .t_initial = 0.008, // initial temperature
    .mu_t = 1.003,
    .t_min = 2.0e-6 // minimum temperature
  };

  gsl_siman_solve (
    rng,
    &xp0,
    ocaml_siman_energy,
    ocaml_siman_step,
    ocaml_siman_distance,
    print_cfg, // print_position
    ocaml_siman_copy, // copy func
    ocaml_siman_construct, // copy cons
    ocaml_siman_destroy, // destructor
    0, // element size
    params
  );

  gsl_rng_free (rng);
  printf ("[ocaml_siman_solve] result: %f\n", Double_val (Field (xp0.state, 0)));
  fflush (stdout);
  CAMLreturn (xp0.state);
}
*/