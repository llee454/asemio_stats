open Core_kernel

let sum ~f = Array.fold ~init:0.0 ~f:(fun sum x -> sum +. f x)

let sumi ~f = Array.foldi ~init:0.0 ~f:(fun i sum x -> sum +. f i x)

let%expect_test "sum_1" =
  printf "%.1f" (sum ~f:Fn.id [| 1.0; 3.5; -2.5; 8.2 |]);
  [%expect {|10.2|}]

external pow_int : float -> int -> float = "ocaml_gsl_pow_int"

let%expect_test "pow_int_1" =
  printf "%.2f" (pow_int 1.1 2);
  [%expect {|1.21|}]

let%expect_test "pow_int_2" =
  printf "%.4f" (pow_int 3.1415 3);
  [%expect {|31.0035|}]

(** Returns x^y *)
let expt (x : float) (y : float) = exp (y *. log x)

let%expect_test "expt_1" =
  printf "%.2f" (expt 3.0 7.4);
  [%expect {|3393.89|}]

external mean : float array -> float = "ocaml_mean"

let%expect_test "mean_1" =
  printf "%.4f" (mean [| 3.1; 2.7; -1.5; 0.5; -3.12 |]);
  [%expect {|0.3360|}]

external fact : int -> float = "ocaml_gsl_sf_fact"

let%expect_test "fact_1" =
  printf "%.1f" (fact 3);
  [%expect {|6.0|}]

let%expect_test "fact_2" =
  printf "%.1f" (fact 2);
  [%expect {|2.0|}]

let%expect_test "fact_3" =
  printf "%.1f" (fact 5);
  [%expect {|120.0|}]

external cdf_gaussian_p : x:float -> std:float -> float = "ocaml_gsl_cdf_gaussian_P"

let%expect_test "cdf_gaussian_p_1" =
  printf "%.4f" (cdf_gaussian_p ~x:0.5 ~std:1.0);
  [%expect {|0.6915|}]

let%expect_test "cdf_gaussian_p_2" =
  printf "%.4f" (cdf_gaussian_p ~x:0.89 ~std:1.0);
  [%expect {|0.8133|}]

let%expect_test "cdf_gaussian_p_3" =
  printf "%.4f" (cdf_gaussian_p ~x:2.14 ~std:3.14);
  [%expect {|0.7522|}]

external cdf_gaussian_q : x:float -> std:float -> float = "ocaml_gsl_cdf_gaussian_Q"

let%expect_test "cdf_gaussian_q_1" =
  printf "%.4f" (1.0 -. cdf_gaussian_q ~x:0.5 ~std:1.0);
  [%expect {|0.6915|}]

let%expect_test "cdf_gaussian_q_2" =
  printf "%.4f" (1.0 -. cdf_gaussian_q ~x:0.89 ~std:1.0);
  [%expect {|0.8133|}]

let%expect_test "cdf_gaussian_q_3" =
  printf "%.4f" (1.0 -. cdf_gaussian_q ~x:2.14 ~std:3.14);
  [%expect {|0.7522|}]

external cdf_chi_squared_p : x:float -> ndf:float -> float = "ocaml_gsl_cdf_chisq_P"

let%expect_test "cdf_chi_squared_p_1" =
  printf "%.4f" (cdf_chi_squared_p ~x:0.0158 ~ndf:1.0);
  [%expect {|0.1000|}]

let%expect_test "cdf_chi_squared_p_2" =
  printf "%.4f" (cdf_chi_squared_p ~x:0.00393 ~ndf:1.0);
  [%expect {|0.0500|}]

let%expect_test "cdf_chi_squared_p_3" =
  printf "%.4f" (cdf_chi_squared_p ~x:0.455 ~ndf:1.0);
  [%expect {|0.5000|}]

let%expect_test "cdf_chi_squared_p_4" =
  printf "%.4f" (cdf_chi_squared_p ~x:3.94 ~ndf:10.0);
  [%expect {|0.0500|}]

external cdf_gaussian_pinv : x:float -> sigma:float -> float = "ocaml_gsl_cdf_gaussian_Pinv"

let%expect_test "cdf_gaussian_pinv_1" =
  printf "%.4f" (cdf_gaussian_pinv ~x:0.5 ~sigma:1.0);
  [%expect {|0.0000|}]

let%expect_test "cdf_gaussian_pinv_2" =
  printf "%.4f" (cdf_gaussian_pinv ~x:0.05 ~sigma:1.0);
  [%expect {|-1.6449|}]

let%expect_test "cdf_gaussian_pinv_2" =
  printf "%.4f" (cdf_gaussian_pinv ~x:0.95 ~sigma:2.0);
  [%expect {|3.2897|}]

external pdf_binomial : k:int -> p:float -> n:int -> float = "ocaml_gsl_ran_binomial_pdf"

let%expect_test "pdf_binomial_1" =
  printf "%.4f" (pdf_binomial ~k:5 ~p:0.5 ~n:10);
  [%expect {|0.2461|}]

module Linear_fit = struct
  type t = {
    c0: float;
    c1: float;
  }

  external ocaml_gsl_fit_linear : float array -> float array -> t = "ocaml_gsl_fit_linear"

  let f = ocaml_gsl_fit_linear
end

module Erf = struct
  external f : float -> float = "ocaml_gsl_sf_erf_Z"

  external q : float -> float = "ocaml_gsl_sf_erf_Q"
end

module Integrate = struct
  type t = {
    out: float;
    err: float;
    neval: int64;
  }

  external ocaml_integrate : f:(float -> float) -> lower:float -> upper:float -> t = "ocaml_integrate"

  let f = ocaml_integrate

  let%expect_test "Integrate.f_1" =
    let { out; _ } =
      f ~lower:(-1.0) ~upper:1.0 ~f:Float.((fun x -> exp ~-(square x /. 2.0) /. (2.0 *. sqrt_pi)))
    in
    printf "%.4f" out;
    [%expect {|0.4827|}]

  external integration_qag : f:(float -> float) -> lower:float -> upper:float -> t
    = "ocaml_integration_qag"

  let g = integration_qag

  let%expect_test "Integrate.g_1" =
    let { out; _ } =
      g ~lower:(-1.0) ~upper:1.0 ~f:Float.((fun x -> exp ~-(square x /. 2.0) /. (2.0 *. sqrt_pi)))
    in
    printf "%.4f" out;
    [%expect {|0.4827|}]

  external integration_qagi : f:(float -> float) -> t = "ocaml_integration_qagi"

  let h = integration_qagi

  let%expect_test "Integrate.h_1" =
    let { out; _ } = h ~f:Float.((fun x -> 2.0 ** ~-(square x))) in
    printf "%.8f" out;
    [%expect {|2.12893404|}]
end

module Nonlinear_fit = struct
  (*
    represents the parameter values used by a function to process an
    argument value
  *)
  type t = {
    ks: float array;
    (* parameter values *)
    x: float; (* argument value *)
  }

  (*
    Accepts a function, f, that has parameters ks and accepts a float
    and returns a float; and a set of initial parameter values,
    ks_init; and a set of data values, xs; and returns the set of
    parameter values that minimize the error between f_ks (x)
  *)
  external f : f:(t -> float) -> ks_init:float array -> xs:float array -> ys:float array -> float array
    = "ocaml_gsl_fit_nlinear"
end

module Stats_tests = struct
  exception N_too_small

  (*
    Computes the Jarque-Bera statistic for the given dataset.

    Note: for n >= 2,000 the distribution of the JB statistic
    approximates a X^2 (Chi-squared) distribution with 2 degrees of
    freedom.

    Warning: this function will throw the N_too_small exception if the
    dataset is smaller than 2,000.

    References:
    1. "The Shapiro-Wilk And Related Tests for Normality",
       [https://math.mit.edu/~rmd/465/shapiro.pdf]
  *)
  let jarque_bera ~xs =
    let n =
      let n = Array.length xs in
      if n < 2000 then raise N_too_small;
      float n
    in
    let u = mean xs in
    (* compute the biased sample deviation *)
    let s = sqrt (sum ~f:(fun x -> Float.square (x -. u)) xs /. n) in
    let skew = sum ~f:(fun x -> pow_int ((x -. u) /. s) 3) xs /. n in
    let kurt = (sum ~f:(fun x -> pow_int ((x -. u) /. s) 4) xs /. n) -. 3.0 in
    n /. 6.0 *. (Float.square skew +. (Float.square kurt /. 4.0))

  (**
    Computes an approximation of the r-th largest value in a random
    sample of n values drawn from a Gaussian distribution having
    the given mean and standard deviation sample.

    Note: this function uses an approximation proposed by Blom (1958).

    Warning: for small n, this approximation is rather rough.

    References:
    1. [Reference Table](https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2Faoms%2F1177728266)
    2. [Algorithm AS 177: Expected Normal Order Statistics (Exact and Approximate)](https://www.jstor.org/stable/2347982 )
  *)
  let order_stat_approx ~r ~n ~mean ~std : float =
    let alpha = 0.375 in
    mean
    +. (cdf_gaussian_pinv ~x:((float r -. alpha) /. (float n -. (2.0 *. alpha) +. 1.0)) ~sigma:1.0 *. std)

  let%expect_test "order_stat_1" =
    printf "%.1f" @@ order_stat_approx ~r:20 ~n:20 ~mean:0.0 ~std:1.0;
    [%expect {|1.9|}]

  (**
    Computes an approximation of mean value of the r-th largest
    values in a random sample of n values drawn from a Gaussian
    distribution with mean 0 and standard deviation 1.

    Note: this function uses a numerical approximation to the exact
    formula as given in [1].

    k * 'integrate (
      x/(1 - cdf (x)) * B_pdf (n, r, 1 - cdf (x)) * pdf (x),
      x);

    References:
    1. [Algorithm AS 177: Expected Normal Order Statistics (Exact and Approximate)](https://www.jstor.org/stable/2347982 )
    2. D. Teichroew, et. al. "Tables of Expected Values of Order
       Statistics and Products of Order Statistics for Samples of
       Size Twenty and Less from the Normal Distribution".
       Numerical Analysis Research, University of California,
       Los Angeles. 
       https://projecteuclid.org/journalArticle/Download?urlId=10.1214%2Faoms%2F1177728266
    3. [Exptected Values of Normal Order Statistics](http://faculty.washington.edu/heagerty/Books/Biostatistics/TABLES/NormalOrder.pdf)
  *)
  let order_stat ~r ~n : float =
    let lower, upper = -6.0, 6.0 in
    (* NOTE: the tests fail when these bounds are expanded. *)
    let Integrate.{ out; _ } =
      Integrate.g ~lower ~upper ~f:(fun x ->
          let p = 1.0 -. Erf.q x in
          float r *. x /. p *. pdf_binomial ~k:r ~p ~n *. Erf.f x)
    in
    out

  let%expect_test "order_stat_1" =
    printf "%.4f" (order_stat ~r:8 ~n:10);
    [%expect {|0.6561|}]

  let%expect_test "order_stat_2" =
    printf "%.4f" (order_stat ~r:2 ~n:4);
    [%expect {|-0.2970|}]

  (*
    Computes the Shapiro Francia statistic which parameterizes
    the likelihood that the given sample is drawn from a Gaussian
    distribution with mean 0 and standard deviation 1.
  *)
  let shapiro_francia_stat (xs : float array) : float =
    let n = Array.length xs in
    let m = mean xs in
    Array.sort ~len:n xs ~compare:[%compare: float];
    let os = Array.init n ~f:(fun r -> order_stat ~r ~n) in
    sumi xs ~f:(fun r y -> (y -. m) *. os.(r))
    /. Float.sqrt
         (sumi xs ~f:(fun _ y -> Float.square (y -. m)) *. sumi xs ~f:(fun r _ -> Float.square os.(r)))

  (*
  let%expect_test "shapiro_francia_stat_1" =
    printf "%.4f"
      (shapiro_francia_stat
         [|
           0.664220;
           -0.631950;
           -0.448290;
           0.184850;
           -1.40500;
           0.896160;
           -0.598050;
           -0.425810;
           0.504560;
           0.732380;
           1.91910;
           -0.0628270;
           0.451500;
           -0.581380;
           -1.07650;
           -0.245060;
           0.204370;
           -0.646910;
           -0.007770;
           -1.47800;
           -0.573960;
           0.448420;
           -1.25420;
           0.220640;
           -1.18590;
           -1.14360;
           -0.890480;
           -0.90406;
           1.24900;
           -0.875340;
           -0.0253890;
           -0.0810450;
           1.87970;
           0.930510;
           0.865260;
           -0.618640;
           0.110770;
           -2.00530;
           0.328060;
           0.593620;
         |]);
    [%expect {||}]
*)

  (*
    Returns the normalized Shapiro Francia Statistics.

    Note: According to P. Royston "A Toolkit for Testing
    for Non-Normality in Complete and Censored Samples", this
    value is normally distributed with mean:

      mean: ln(ln(n)) - ln(n)
      std: ln(ln(n)) + 2/ln(n)
    
    where n denotes the sample size.

    This transformation was originally reported by Roysten (1992a)
    for n values from 5 to 5000.

    Warning: this transformation is only valid for n values ranging
    from 5 to 5000.

    References:
    1. ["A Toolkit for Testing for Non-Normality in Complete and Censored Samples"]( https://www.jstor.org/stable/2348109)
  *)

  let shapiro_francia_test (xs : float array) : float =
    let n = float (Array.length xs) in
    let mean = Float.log (Float.log n) -. Float.log n in
    let std = Float.log (Float.log n) +. (2.0 /. Float.log n) in
    let w = Float.log (1.0 -. shapiro_francia_stat xs) in
    let p = cdf_gaussian_p ~x:(w -. mean) ~std in
    let q = 1.0 -. p in
    if Float.(p >= 0.5) then p -. q else q -. p

  (*
  let%expect_test "shapiro_francia_test" =
    printf "%.4f"
      (shapiro_francia_test
         [|
           0.664220;
           -0.631950;
           -0.448290;
           0.184850;
           -1.40500;
           0.896160;
           -0.598050;
           -0.425810;
           0.504560;
           0.732380;
           1.91910;
           -0.0628270;
           0.451500;
           -0.581380;
           -1.07650;
           -0.245060;
           0.204370;
           -0.646910;
           -0.007770;
           -1.47800;
           -0.573960;
           0.448420;
           -1.25420;
           0.220640;
           -1.18590;
           -1.14360;
           -0.890480;
           -0.90406;
           1.24900;
           -0.875340;
           -0.0253890;
           -0.0810450;
           1.87970;
           0.930510;
           0.865260;
           -0.618640;
           0.110770;
           -2.00530;
           0.328060;
           0.593620;
         |]);
    [%expect {||}]
*)

  (*
  let%expect_test "shapiro_francia_test_2" =
    let () = Random.init 5 in
    let xs = Array.init 30 ~f:(fun _ -> Random.float 3.0) in
    printf "%.4f" (shapiro_francia_test xs);
    [%expect {||}]
*)
end

module Simulated_annealing = struct
  type 'a t = {
    cons: 'a t ref -> 'a t ref;
    energy: 'a ref -> float;
    step: 'a t ref -> max_dist:float -> unit;
    dist: 'a ref -> 'a ref -> float;
    state: 'a ref;
    copy: source:'a t ref -> dest:'a t ref -> unit;
  }
  [@@deriving sexp]

  external simulated_annealing : 'a t ref -> 'a ref = "ocaml_siman_solve"

  let f = simulated_annealing
end
