//
// cli.rs
// Defines a simulation with no visual cues.
// Is not accessed if no command-line arguments are given - the code is diverted to gui.rs instead in such a case.
//

use std::env;
use sim::simulation;
use time;

pub fn parse_args(){
    let mut valid = true;
    let will_print: bool;
    let mut clargs = vec![];
    let start = time::PreciseTime::now();
    for argument in env::args() {
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
        let results = simulation(clargs[0],
                    clargs[1],
                    clargs[2],
                    clargs[3],
                    clargs[4],
                    clargs[5],
                    clargs[6],
                    clargs[7],
                    clargs[8] as i32,
                    will_print);
        let mut fails = 0;
        let mut avg_successes = 0.;
        for result in results {
            match result {
                Some(hits) => {
                    avg_successes += hits as f64;
                },
                None => {
                    fails += 1;
                }
            }
        }
        avg_successes = avg_successes/(clargs[8]-fails as f64);
        println!("{} rays passed through the fiber.",(clargs[8] as i32)-fails);
        println!("{} rays returned to fiber entrance or got lost", fails);
        println!("{} is the average number of hits among successful passes", avg_successes);
        println!("{:?} execution time", (start.to(time::PreciseTime::now())) );
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
