//
// cli.rs
// Defines a simulation with no visual cues.
// Is not accessed if no command-line arguments are given - the code is diverted to gui.rs instead in such a case.
//

use std::env;

pub fn parse_args(){
    //
    for argument in env::args() {
        // Here we want to actually parse the CL arguments and not print them
        println!("{}", argument);
        //
    }
}
