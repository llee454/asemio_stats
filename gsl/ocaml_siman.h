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
  type 'a t = {
    copy: 'a ref -> 'a ref;
    energy: 'a ref -> float;
    step: 'a ref -> float -> unit;
    dist: 'a ref -> 'a ref -> float;
    state: 'a ref
  }
*/

void ocaml_siman_print_context (value xp) {
  CAMLparam1 (xp);
  printf ("[ocaml_siman_print_context]\n");
  printf (
    "{ cons: %lu; energy: %lu; step: %lu; dist: %lu; state: %lu; copy: %lu }\n",
    Field (Field (xp, 0), 0),
    Field (Field (xp, 0), 1),
    Field (Field (xp, 0), 2),
    Field (Field (xp, 0), 3),
    Field (Field (xp, 0), 4),
    // Double_val (Field (Field (Field (xp, 0), 4), 0)),
    Field (Field (xp, 0), 5)
  );
  fflush (stdout);
  CAMLreturn0;
}

double ocaml_siman_energy (value xp) {
  CAMLparam1 (xp);
  printf ("[ocaml_siman_energy]\n");
  fflush (stdout);
  printf ("[ocaml_siman_energy] xp: %lu\n", xp);
  fflush (stdout);
  ocaml_siman_print_context (xp);
  double result = Double_val (caml_callback (Field (Field (xp, 0), 1), Field (Field (xp, 0), 4)));
  ocaml_siman_print_context (xp);
  CAMLreturnT (double, result);
}

void ocaml_siman_step (const gsl_rng* rng, value xp, double step_size) {
  CAMLparam1 (xp);
  printf ("[ocaml_siman_step]\n");
  fflush (stdout);
  printf ("[ocaml_siman_step] xp: %lu\n", xp);
  fflush (stdout);
  ocaml_siman_print_context (xp);
//   double step   = (2 * gsl_rng_uniform (rng) - 1) * step_size;
  caml_callback2 (Field (Field (xp, 0), 2), xp, caml_copy_double (step_size));
  printf ("[ocaml_siman_step] after step\n");
  ocaml_siman_print_context (xp);
  CAMLreturn0;
}

double ocaml_siman_distance (value xp, value yp) {
  CAMLparam2 (xp, yp);
  printf ("[ocaml_siman_distance]\n");
  fflush (stdout);
  printf ("[ocaml_siman_distance] xp: %lu yp: %lu\n", xp, yp);
  fflush (stdout);
  CAMLreturnT (double, Double_val (caml_callback2 (Field (Field (xp, 0), 3), Field (Field (xp, 0), 4), Field (Field (yp, 0), 4))));
}

void ocaml_siman_copy (value source, value dest) {
  CAMLparam2 (source, dest);
  printf ("[ocaml_siman_copy]\n");
  fflush (stdout);
  void* orig_source = (void*) source;
  void* orig_dest   = (void*) dest;
  printf ("[ocaml_siman_copy] initial source: %lu initial dest: %lu\n", source, dest);
  fflush (stdout);
  ocaml_siman_print_context (source);
  ocaml_siman_print_context (dest);
  caml_callback2 (Field (Field (source, 0), 5), source, dest);
  if (orig_source != (void*) source || orig_dest != (void*) dest) {
    printf ("[ERROR] source or dest reference pointer changed after/during call to copy callback.\n");
    fflush (stdout);
  }
  printf ("[ocaml_siman_copy] source: %lu dest: %lu source fn: %lu\n", source, dest, Field (Field (source, 0), 5));
  ocaml_siman_print_context (source);
  ocaml_siman_print_context (dest);
  CAMLreturn0;
}

value ocaml_siman_construct (value xp) {
  CAMLparam1 (xp);
  printf ("[ocaml_siman_construct]\n");
  fflush (stdout);
  CAMLlocal1 (result);
  result = caml_callback (Field (Field (xp, 0), 0), xp);
  printf ("[ocaml_siman_construct] xp: %lu result: %lu value: %0.4f\n", xp, result, Double_val (Field (Field (Field (result, 0), 4), 0)));
  ocaml_siman_print_context (xp);
  ocaml_siman_print_context (result);
  fflush (stdout);
  CAMLreturn (result);
}

void ocaml_siman_destroy (value xp) {
  CAMLparam0 ();
  printf ("[ocaml_siman_destroy]\n");
  fflush (stdout);
  CAMLreturn0;
}

void ocaml_siman_print (value xp) {
  CAMLparam0 ();
  printf (" x = %0.4f", Double_val (Field (Field (Field (xp, 0), 4), 0)));
  fflush (stdout);
  CAMLreturn0;
}

CAMLprim value ocaml_siman_solve (value context) {
  CAMLparam1 (context);
  printf ("[ocaml_siman_solve] context: %lu\n", context);
  fflush (stdout);

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

  printf ("[ocaml_siman_solve] done\n");
  fflush (stdout);

  gsl_rng_free (rng);
  CAMLreturn (Field (Field (context, 0), 4));
}