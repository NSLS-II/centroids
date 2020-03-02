Second Pass Through Photons
===========================

Once potential photon events have been located, the centroiding algorithm
processes the list of potential events. If an event is determined to be a
photon, then the resultant parameters are calculated and added to a list of
detected photons. The flow of this algorithm is shown in
:ref:`this flowchart<second_pass_digraph>`.

The algorithm starts by assessing what the potential photon's overlap with
other events is. This is done using the overlap image calculated in the
firstpass. Here if the sum of the :math:`n \times n` cluster is calculated.
If the sum is equal to the cluster size, then the photon event is isolated.
If the value is greater than that of the cluster size, then the photon event
overlaps another event. If this value is greater than the value specified in
the parameters (:cpp:var:`overlap_max`), then the event is rejected.

Next, the pixel intensities in the :math:`n \times n` box are sorted, highest
first via a `bubble sort <https://en.wikipedia.org/wiki/Bubble_sort>`_ . The
background value is calculated by averaging the last
:cpp:var:`pixel_photon_bgnd` values in the sorted list. The integrated
intensity of the photon is then calculated by summing up the first
:cpp:var:`pixel_photon_num` values is the sorted list after subtracting the
average background value just calculated. The photon event is then filtered
using the integrated intensity by the :cpp:var:`sum_min` and
:cpp:var:`sum_max`.

After calculating the integrated intensity, the *center of mass* of the
photon event is calculated using the equations:

.. math::

   C_x = \frac{\sum_i I_i x_i}{\sum_i I_i}
   \quad\text{and}\quad
   C_y = \frac{\sum_i I_i y_i}{\sum_i I_i}

where :math:`I_x` and :math:`I_y` are the intensity of the background
corrected pixel intensities, :math:`x` and :math:`y` are the fractional
pixel coordinates and the sum index :math:`i` is over all the pixels in the
cluster.

The pixel cluster is then fitted using a least squares algorithm in both 2D
and 1D. The pixel cluster is fitted to the *error function* being the
integral of a *gaussian function*. In 2D this function is:

.. math::
   \DeclareMathOperator\erf{erf}

   I = B + I_0 \left(
   \erf \frac{\left( x - \frac{1}{2} - x_0 \right)}{\sqrt{2} \sigma} -
   \erf \frac{\left( x + \frac{1}{2} - x_0 \right)}{\sqrt{2} \sigma}
   \right)\cdot\left(
   \erf \frac{\left( y - \frac{1}{2} - y_0 \right)}{\sqrt{2} \sigma} -
   \erf \frac{\left( y + \frac{1}{2} - y_0 \right)}{\sqrt{2} \sigma}
   \right)


With the 1D version being:

.. math::
   \DeclareMathOperator\erf{erf}

   I = B + I_0 \left(
   \erf \frac{\left( x - \frac{1}{2} - x_0 \right)}{\sqrt{2} \sigma} -
   \erf \frac{\left( x + \frac{1}{2} - x_0 \right)}{\sqrt{2} \sigma}
   \right)


where :math:`I` is the intensity at pixel :math:`x` (or :math:`x,y` in 2D),
:math:`B` is the background and :math:`x_0` and :math:`y_0` are the center
coordinates of the gaussian.

The results of the fit are stored in the output table of photon parameters.

.. _second_pass_digraph:

.. digraph:: second_pass
   :caption: Second pass through events found in first pass

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