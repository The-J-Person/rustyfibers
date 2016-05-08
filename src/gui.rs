//
// gui.rs
// Defines and manages the GUI components of the program.
// Is not accessed if any command-line arguments are given - the code is diverted to cli.rs instead in such a case.
//

extern crate gtk;

use gtk::prelude::*;

fn clicked_exit(button: &gtk::Button){
    //End all threads...?
    gtk::main_quit();
    Inhibit(false);
}

pub fn start_gui() {
    if gtk::init().is_err() {
        println!("Failed to initialize GTK.");
        return;
    }

    let window = gtk::Window::new(gtk::WindowType::Toplevel);

    window.set_title("Rusty Ray Tracer");
    window.set_border_width(10);
    window.set_position(gtk::WindowPosition::Center);
    //window.set_default_size(350, 70);

    window.connect_delete_event(|_, _| {
        gtk::main_quit();
        Inhibit(false)
    });

    //

    // Button Definitions
    let button_start = gtk::Button::new_with_label("Start");
    let button_pause = gtk::Button::new_with_label("Pause");
    let button_help = gtk::Button::new_with_label("Help");
    let button_exit = gtk::Button::new_with_label("Exit");

    button_exit.connect_clicked(clicked_exit);

    // Layout Element Definitions
    let frame_properties = gtk::Frame::new(Some("Properties"));
    frame_properties.set_border_width(10);
    let frame_main = gtk::Frame::new(Some("Main"));
    frame_main.set_border_width(10);
    let frame_results = gtk::Frame::new(Some("Results"));
    frame_results.set_border_width(10);
        // False Video spaceholder
    let frame_video = gtk::Frame::new(Some("Video"));
    frame_video.set_border_width(10);

    let h_box = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let v_box_left = gtk::Box::new(gtk::Orientation::Vertical, 10);
    let v_box_right = gtk::Box::new(gtk::Orientation::Vertical, 10);
    let h_box_m = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let v_box_m = gtk::Box::new(gtk::Orientation::Vertical, 10);
    let h_box_m1 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_m2 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p01 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p02 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p03 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p04 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p05 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p06 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p07 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p08 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p09 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p10 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p11 = gtk::Box::new(gtk::Orientation::Horizontal, 10);
    let h_box_p12 = gtk::Box::new(gtk::Orientation::Horizontal, 10);

    let label_n1 = gtk::Label::new(Some("n₁"));
    let entry_n1 = gtk::Entry::new();
    entry_n1.set_text("1.335");
    let label_n2 = gtk::Label::new(Some("n₂"));
    let entry_n2 = gtk::Entry::new();
    entry_n2.set_text("1.457");
    let label_radius = gtk::Label::new(Some("Radius[mm]:"));
    let scale_radius = gtk::SpinButton::new_with_range(0.05, 5., 0.05);
    let label_length = gtk::Label::new(Some("Fiber Length[cm]:"));
    let entry_length = gtk::Entry::new();
    entry_length.set_text("10");
    let label_length_untapped = gtk::Label::new(Some("Untapped Length[cm]:"));
    let entry_length_untapped = gtk::Entry::new();
    entry_length_untapped.set_text("5.0");
    let label_radius_untapped = gtk::Label::new(Some("Untapped Radius[mm]:"));
    let scale_radius_untapped = gtk::SpinButton::new_with_range(0.05, 5., 0.05);
    //let check_change_angle = gtk::CheckButton::new_with_label("Change Angle");
    let label_angle = gtk::Label::new(Some("Angle[c]:"));
    let scale_angle = gtk::SpinButton::new_with_range(0., 90., 0.05);
    let check_diamond = gtk::CheckButton::new_with_label("Diamond");
    let label_histories = gtk::Label::new(Some("Histories:"));
    let entry_histories = gtk::Entry::new();
    let check_fast_histories = gtk::CheckButton::new_with_label("Fast");
    let check_show_fiber = gtk::CheckButton::new_with_label("Show Fiber");
    let check_show_animation = gtk::CheckButton::new_with_label("Show Animation");

    let tree_prop = gtk::TreeView::new();
    let treecol_p1 = gtk::TreeViewColumn::new();
    treecol_p1.set_title("ID");
    tree_prop.append_column(&treecol_p1);
    let treecol_p2 = gtk::TreeViewColumn::new();
    treecol_p2.set_title("Status");
    tree_prop.append_column(&treecol_p2);
    let treecol_p3 = gtk::TreeViewColumn::new();
    treecol_p3.set_title("Progress");
    tree_prop.append_column(&treecol_p3);
    let treecol_p4 = gtk::TreeViewColumn::new();
    treecol_p4.set_title("Hit Progress");
    tree_prop.append_column(&treecol_p4);
    let treecol_p5 = gtk::TreeViewColumn::new();
    treecol_p5.set_title("No Hits");
    tree_prop.append_column(&treecol_p5);

    let tree_results = gtk::TreeView::new();
    let treecol_r1 = gtk::TreeViewColumn::new();
    treecol_r1.set_title("Passed");
    tree_results.append_column(&treecol_r1);
    let treecol_r2 = gtk::TreeViewColumn::new();
    treecol_r2.set_title("Failed");
    tree_results.append_column(&treecol_r2);
    let treecol_r3 = gtk::TreeViewColumn::new();
    treecol_r3.set_title("Hits");
    tree_results.append_column(&treecol_r3);
    let treecol_r4 = gtk::TreeViewColumn::new();
    treecol_r4.set_title("Ratio %");
    tree_results.append_column(&treecol_r4);

    let label_midcenter = gtk::Label::new(Some("Midcenter Type:"));
    let combo_midcenter = gtk::ComboBoxText::new();
    combo_midcenter.append_text("Radial");
    combo_midcenter.append_text("Flattened");
    combo_midcenter.set_active(0);

    // Layout element attachment
    window.add(&h_box);
    h_box.add(&frame_properties);
    h_box.add(&v_box_right);
    frame_properties.add(&v_box_left);
    v_box_right.add(&h_box_m);
    v_box_right.add(&frame_video);
    h_box_m.add(&frame_main);
    frame_main.add(&v_box_m);
    h_box_m.add(&frame_results);
    frame_results.add(&tree_results);
    tree_results.columns_autosize();
    v_box_m.add(&h_box_m1);
    v_box_m.add(&h_box_m2);
    h_box_m2.add(&button_start);
    h_box_m2.add(&button_pause);
    h_box_m2.add(&button_help);
    h_box_m2.add(&button_exit);
    v_box_left.add(&h_box_p01);
    h_box_p01.add(&label_n1);
    h_box_p01.add(&entry_n1);
    h_box_p01.add(&label_n2);
    h_box_p01.add(&entry_n2);
    v_box_left.add(&h_box_p02);
    h_box_p02.add(&label_radius);
    h_box_p02.add(&scale_radius);
    v_box_left.add(&h_box_p03);
    h_box_p03.add(&label_length);
    h_box_p03.add(&entry_length);
    v_box_left.add(&h_box_p04);
    h_box_p04.add(&label_length_untapped);
    h_box_p04.add(&entry_length_untapped);
    v_box_left.add(&h_box_p05);
    h_box_p05.add(&label_radius_untapped);
    h_box_p05.add(&scale_radius_untapped);
    v_box_left.add(&h_box_p06);
    h_box_p06.add(&label_angle);
    h_box_p06.add(&scale_angle);
    v_box_left.add(&h_box_p07);
    h_box_p07.add(&check_diamond);
    v_box_left.add(&h_box_p08);
    h_box_p08.add(&label_histories);
    h_box_p08.add(&entry_histories);
    h_box_p08.add(&check_fast_histories);
    v_box_left.add(&h_box_p09);
    h_box_p09.add(&check_show_fiber);
    v_box_left.add(&h_box_p10);
    h_box_p10.add(&check_show_animation);
    v_box_left.add(&h_box_p11);
    h_box_p11.add(&tree_prop);
    tree_prop.columns_autosize();
    v_box_left.add(&h_box_p12);
    h_box_p12.add(&label_midcenter);
    h_box_p12.add(&combo_midcenter);

    window.show_all();
    gtk::main();
}
