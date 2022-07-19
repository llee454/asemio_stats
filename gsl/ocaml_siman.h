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

#include <gsl_siman.h> // simulated annealing
#include <gsl_rng.h> // random number generator

/*
  type 'a context = {
    copy: 'a ref -> 'a ref,
    energy: 'a ref -> float,
    step: 'a ref -> float -> unit,
    dist: 'a ref -> 'a ref -> float,
    state: 'a ref
  }
*/

double ocaml_siman_energy (value xp) {
  CAMLparam1 (xp);
  CAMLreturnT (double, Double_val (caml_callback (Field (xp, 1), Field (xp, 4))));
}

void ocaml_siman_step (const gsl_rng* rng, value xp, double step_size) {
  CAMLparam1 (xp);
  double step   = (2 * gsl_rng_uniform (rng) - 1) * step_size;
  caml_callback2 (Field (xp, 2), Field (xp, 4), caml_copy_double (step));
  CAMLreturn0;
}

double ocaml_siman_distance (value xp, value yp) {
  CAMLparam2 (xp, yp);
  CAMLreturnT (double, Double_val (caml_callback2 (Field (xp, 3), Field (xp, 4), Field (yp, 4))));
}

void ocaml_siman_copy (value source, value dest) {
  CAMLparam2 (source, dest);
  CAMLreturn0;
}

value ocaml_siman_construct (value xp) {
  CAMLparam1 (xp);
  CAMLlocal1 (result);
  result = caml_alloc(5, 0);
  Store_field (result, 0, Field (xp, 0));
  Store_field (result, 1, Field (xp, 1));
  Store_field (result, 2, Field (xp, 2));
  Store_field (result, 3, Field (xp, 3));
  Store_field (result, 4, caml_callback (Field (xp, 1), Field (xp, 4)));
  CAMLreturn (result);
}

void ocaml_siman_destroy (value xp) {
  CAMLparam0 ();
  CAMLreturn0;
}

void ocaml_siman_print (value xp) {
  CAMLparam0 ();
  CAMLreturn0;
}

CAMLprim value ocaml_siman_solve (value context) {
  CAMLparam1 (context);

  gsl_rng_env_setup ();
  gsl_rng* rng = gsl_rng_alloc (gsl_rng_default);

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
    (void*) context,
    (gsl_siman_Efunc_t) ocaml_siman_energy,
    (gsl_siman_step_t) ocaml_siman_step,
    (gsl_siman_metric_t) ocaml_siman_distance,
    (gsl_siman_print_t) ocaml_siman_print,
    (gsl_siman_copy_t) ocaml_siman_copy,
    (gsl_siman_copy_construct_t) ocaml_siman_construct,
    (gsl_siman_destroy_t) ocaml_siman_destroy,
    0, // element size
    params
  );

  gsl_rng_free (rng);
  CAMLreturn (Field (context, 4));
}