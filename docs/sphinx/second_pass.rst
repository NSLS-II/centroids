Second Pass Through Photons
===========================

Second pass through image

.. digraph:: second_pass

    label="Initial Pass Through Image"

    node[shape="box", style=rounded]
       start; end;
    node[shape="diamond", style=""]
       is_last_event[label="is the last event?"]
       is_overlap_ok[label="is the event overlapping\nmore than specified?"]
       is_int_ok[label="is the integral\nwithin bounds?"]
    node[shape="parallelogram", style=""]
       calc_overlap[label="calculate overlap"]
       sort_pixels[label="sort pixels by intensity\nin event nxn box"]
       calc_bgnd[label="calculate background"]
       calc_int[label="calculate integral"]
       calc_com[label="calculate com"]
       calc_2d[label="calculate 2D fit"]
       calc_1d[label="calculate 1D fits"]
       store[label="store values in photon table"]
    node[shape="box", style=""]
       set_start_event[label="move to first event in list"]
       move_to_next_event[label="move to next event"]

    start -> set_start_event
    set_start_event -> calc_overlap
    calc_overlap -> is_overlap_ok
    is_overlap_ok -> is_last_event[label="> overlap_max"]
    is_overlap_ok -> sort_pixels[label="< overlap_max"]
    sort_pixels -> calc_bgnd
    calc_bgnd -> calc_int
    calc_int -> is_int_ok
    is_int_ok -> is_last_event[label="out of bounds"]
    is_int_ok -> calc_com[label="sum_min < int < sum_max"]
    calc_com -> calc_2d
    calc_2d -> calc_1d
    calc_1d -> store
    store -> is_last_event

    is_last_event -> end[label="yes"]
    is_last_event -> move_to_next_event[label="no"]
    move_to_next_event -> calc_overlap