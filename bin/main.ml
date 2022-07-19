open! Core_kernel
open! Lwt.Syntax
open! Lwt.Infix
open! Gsl

let main =
  Lwt_main.run @@
    let Linear_fit.{c0; c1} = Linear_fit.f [|0.0; 1.0; 2.0|] [|5.0; 10.0; 20.0|] in 
    let* () = Lwt_io.printlf "linear regression: %f + %f x" c0 c1 in
    let Integrate.{ out; err; _ } = Integrate.f ~f:(fun x -> Float.log1p x) ~lower:0.1 ~upper:0.5 in
    let* () = Lwt_io.printlf "integration: result %f err %f" out err in
    let result =
      Simulated_annealing.(f
        {
        copy = (fun x -> print_endline @@ sprintf "[copy] %f" !x;ref !x);
        energy = (fun x -> (!x +. 1.0) *. (!x +. 1.0));
        step = (fun x dist -> x := !x +. dist);
        dist = (fun x y -> Float.abs (!x -. !y));
        state = (ref 5.0);
        })
    in
    let* () = Lwt_io.printlf "simulated annealing result: %f" !result in
    Lwt.return_unit
