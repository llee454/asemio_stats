#include <stdio.h>

// Prevent OCaml from exporting short macro names.
#define CAML_NAME_SPACE 1

#include <caml/callback.h>
#include <caml/mlvalues.h>
#include <caml/memory.h> // CAMLreturn
#include <caml/alloc.h> // caml_copy

#include <gsl_math.h> // log1p
#include <gsl_sf_erf.h>
#include <gsl_fit.h> // gsl_fit_linear
#include <gsl_integration.h>


CAMLprim value ocaml_log1p (value x) {
  CAMLparam0 ();
  CAMLreturn (caml_copy_double (gsl_log1p (Double_val (x))));
}

CAMLprim value ocaml_gsl_sf_erf_Z (value x) {
  CAMLparam1 (x);
  CAMLreturn (caml_copy_double (gsl_sf_erf_Z (Double_val (x))));
}

CAMLprim value ocaml_gsl_fit_linear (value xs, value ys) {
  CAMLparam2 (xs, ys);
  CAMLlocal1 (result);
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
    .params = &params
  };
  int status = gsl_integration_qng (&F, Double_val (lower), Double_val (upper), 0.1, 0.1, &out, &err, &neval);
  result = caml_alloc (3, 0);
  Store_field (result, 0, caml_copy_double (out));
  Store_field (result, 1, caml_copy_double (err));
  Store_field (result, 2, Val_long (neval));
  CAMLreturn (result);
}