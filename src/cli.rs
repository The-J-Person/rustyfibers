//
// cli.rs
// Defines a simulation with no visual cues.
// Is not accessed if no command-line arguments are given - the code is diverted to gui.rs instead in such a case.
//

use std::env;
use sim::simulation;

// pub fn simulation(out_len: f64, in_len: f64, cone_len: f64, n1: f64, n2: f64, out_rad: f64,
//                     in_rad: f64, angle: f64, thread_amount: i32, output_to_console: bool)

pub fn parse_args(){
    let mut valid = true;
    let will_print: bool;
    let mut clargs = vec![];
    for argument in env::args() {
        // Here we want to actually parse the CL arguments and not print them
        let n = argument.parse::<f64>();
        match n {
            Ok(num) => {
                clargs.push(num);
            },
            Err(_) => {
                if clargs.len()>0 {
                    valid = false;
                }
            },
        }
    }
    if clargs.len()<10 {
        valid = false;
    }
    if valid {
        if clargs[9]==1. {
            will_print = true;
        }
        else {
            will_print = false;
        }
        simulation(clargs[0],
                    clargs[1],
                    clargs[2],
                    clargs[3],
                    clargs[4],
                    clargs[5],
                    clargs[6],
                    clargs[7],
                    clargs[8] as i32,
                    will_print);
    }
    else {
        print_help();
    }
}

pub fn print_help(){
    println!("Rusty Ray Tracer CLI usage instructions:
            \tYou must these arguments in the order in which they appear:\n
            Length of outer fiber
            Length of inner fiber
            Length of connector cones
            n1
            n2
            Radius of outer fiber
            Radius of inner fiber
            Angle bend of inner fiber (in degrees)
            Number of simulations to run
            0 if only the final result should be outputted, 1 for output on individual simulations.");
}
