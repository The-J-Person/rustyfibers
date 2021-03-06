//
// main.rs
// Simply diverts the code to either CLI or GUI, depending on wether command line arguments are given.
//

#![crate_type = "bin"]

extern crate gtk;
extern crate roots;
extern crate rand;
extern crate nalgebra;
extern crate time;

pub mod gui;
pub mod cli;
use std::env;

pub mod math;
pub mod fiber;
pub mod sim;

fn main() {
    if env::args().count() > 1 {
        cli::parse_args();
    }
    else {
        gui::start_gui();
    }
}
