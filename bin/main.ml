open! Core_kernel
open! Lwt.Syntax
open! Lwt.Infix
open! Gsl

let main =
  Lwt_main.run
  @@
  let Linear_fit.{ c0; c1 } = Linear_fit.f [| 0.0; 1.0; 2.0 |] [| 5.0; 10.0; 20.0 |] in
  let* () = Lwt_io.printlf "linear regression: %f + %f x" c0 c1 in
  let Integrate.{ out; err; _ } = Integrate.f ~f:(fun x -> Float.log1p x) ~lower:0.1 ~upper:0.5 in
  let* () = Lwt_io.printlf "integration: result %f err %f" out err in
  let result =
    let open Simulated_annealing in
    let context =
      {
        cons =
          (fun x_ref ->
            let ({ state = x_state_ref; _ } as x) = !x_ref in
            ref { x with state = ref !x_state_ref });
        energy =
          (fun x ->
            print_endline "[energy]";
            let energy = (!x +. 1.0) *. (!x +. 1.0) in
            print_endline @@ sprintf "[energy] x = %f energy = %f" !x energy;
            energy);
        step =
          (fun x_ref ~max_dist ->
            let { state = state_ref; _ } = !x_ref in
            state_ref := !state_ref +. Random.float max_dist);
        dist =
          (fun x y ->
            print_endline @@ sprintf "[dist] x = %f y = %f" !x !y;
            Float.abs (!x -. !y));
        state = ref 5.0;
        copy =
          (fun ~source:{ contents = { state; _ } as source } ~dest ->
            print_endline @@ sprintf "[copy]";
            dest := { source with state = ref !state };
            print_endline
            @@ sprintf
                 !"[copy] source = %{Sexp} dest = %{Sexp} (precopy)"
                 ([%sexp_of: float t] source)
                 ([%sexp_of: float t] !dest));
      }
    in
    f (ref context)
  in
  let* () = Lwt_io.printlf "simulated annealing result: %f" !result in
  Lwt.return_unit
