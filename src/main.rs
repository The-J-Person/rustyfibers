//
// main.rs
// Simply diverts the code to either CLI or GUI, depending on wether command line arguments are given.
//

#![crate_type = "bin"]

extern crate gtk;

pub mod gui;
pub mod cli;
use std::env;

fn main() {
    if env::args().count() > 1 {
        cli::parse_args();
    }
    else {
        gui::start_gui();
    }
}