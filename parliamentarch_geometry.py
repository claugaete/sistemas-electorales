"""
El siguiente código es un fragmento de:
https://github.com/Gouvernathor/parliamentarch
"""

"""
BSD 3-Clause License

Copyright (c) 2024, Gouvernathor

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


import enum
import functools
import math

__all__ = ("get_row_thickness", "get_rows_from_nrows", "get_nrows_from_nseats", "FillingStrategy", "get_seats_centers")

# default angle, in degrees, coming from the rightmost seats through the center to the leftmost seats
_DEFAULT_SPAN_ANGLE = 180

def get_row_thickness(nrows: int) -> float:
    """
    Returns the thickness of a row in the same unit as the coordinates.
    """
    return 1 / (4*nrows - 2)

def get_rows_from_nrows(nrows: int, span_angle: float = _DEFAULT_SPAN_ANGLE) -> list[int]:
    """
    This indicates the maximal number of seats for each row for a given number of rows.
    Returns a list of number of seats per row, from inner to outer.
    The length of the list is nrows.
    span_angle, if provided, is the angle in degrees that the hemicycle, as an annulus arc, covers.
    """
    rv = []

    # thickness of a row (as an annulus) compared to the outer diameter of the hemicycle
    # this is equal to the diameter of a single seat
    rad = get_row_thickness(nrows)
    # if you divide the half-disk of the hemicycle
    # into one half-disk of half the radius
    # and one half-annulus outside it,
    # the innermost row lies on the border between the two,
    # and the outermost row lies entirely inside the half-annulus.
    # So, looking at the line cutting the circle and the annulus in half
    # (which is the bottom border of the diagram),
    # all rows minus one half of the innermost are on the left, same on the right,
    # and the radius of the void at the center is equal to that value again.
    # So, total = 4 * (nrows-.5) = 4*nrows - 2

    radian_span_angle = math.pi*span_angle/180

    for r in range(nrows):
        # row radius : the radius of the circle crossing the center of each seat in the row
        row_arc_radius = .5 + 2*r*rad

        rv.append(int(radian_span_angle*row_arc_radius/(2*rad)))

    return rv

@functools.cache
def _cached_get_rows_from_nrows(nrows: int, span_angle: float = _DEFAULT_SPAN_ANGLE) -> tuple[int, ...]:
    """
    Returns tuples to avoid cache mutation issues.
    """
    return tuple(get_rows_from_nrows(nrows, span_angle))

def get_nrows_from_nseats(nseats: int, span_angle: float = _DEFAULT_SPAN_ANGLE) -> int:
    """
    Returns the minimal number of rows necessary to contain nseats seats.
    """
    i = 1
    while sum(_cached_get_rows_from_nrows(i, span_angle)) < nseats:
        i += 1
    return i

_cached_get_nrows_from_nseats = functools.cache(get_nrows_from_nseats)


class FillingStrategy(enum.StrEnum):
    def __new__(cls, value, doc):
        self = str.__new__(cls, value)
        self._value_ = value
        self.__doc__ = doc
        return self

    DEFAULT = enum.auto(), """
    The seats are distributed among all the rows,
    proportionally to the maximum number of seats each row may contain.
    """

    EMPTY_INNER = enum.auto(), """
    Selects as few outermost rows as necessary, then distributes the seats among them,
    proportionally to the maximum number of seats each row may contain.
    """

    OUTER_PRIORITY = enum.auto(), """
    Fills up the rows as much as possible, starting with the outermost ones.
    """

def get_seats_centers(nseats: int, *,
                      min_nrows: int = 0,
                      filling_strategy: FillingStrategy = FillingStrategy.DEFAULT,
                      span_angle: float = _DEFAULT_SPAN_ANGLE,
                      ) -> dict[tuple[float, float], float]:
    """
    Returns a list of nseats seat centers as (angle, x, y) tuples.
    The canvas is assumed to be of 2 in width and 1 in height, with the y axis pointing up.
    The angle is calculated from the (1., 0.) center of the hemicycle, in radians,
    with 0° for the leftmost seats, 90° for the center and 180° for the rightmost.

    The minimum number of rows required to contain the given number of seats
    will be computed automatically.
    If min_nrows is higher, that will be the number of rows, otherwise the parameter is ignored.
    Passing a higher number of rows will make the diagram sparser.

    seat_radius_factor should be between 0 and 1,
    with seats touching their neighbors in packed rows at seat_radius_factor=1.
    It is only taken into account when placing the seats at the extreme left and right of the hemicycle
    (which are the seats at the bottom of the diagram),
    although the placement of these seats impacts in turn the placement of the other seats.

    span_angle is the angle in degrees from the rightmost seats,
    through the center, to the leftmost seats.
    It defaults to 180° to make a true hemicycle.
    Values above 180° are not supported.
    """
    nrows = max(min_nrows, _cached_get_nrows_from_nseats(nseats, span_angle))
    # thickness of a row in the same unit as the coordinates
    row_thicc = get_row_thickness(nrows)
    span_angle_margin = (1 - span_angle/180)*math.pi /2

    maxed_rows = _cached_get_rows_from_nrows(nrows, span_angle)

    match filling_strategy:
        case FillingStrategy.DEFAULT:
            starting_row = 0
            filling_ratio = nseats/sum(maxed_rows)

        case FillingStrategy.EMPTY_INNER:
            rows = list(maxed_rows)
            while sum(rows[1:]) >= nseats:
                rows.pop(0)
            # here, rows represents the rows which are enough to contain nseats,
            # and their number of seats

            # this row will be the first one to be filled
            # the innermore ones are empty
            starting_row = nrows-len(rows)
            filling_ratio = nseats/sum(rows)
            del rows

        case FillingStrategy.OUTER_PRIORITY:
            rows = list(maxed_rows)
            while sum(rows) > nseats:
                rows.pop(0)
            # here, rows represents the rows which will be fully filled,
            # and their number of seats

            # this row will be the only one to be partially filled
            # the innermore ones are empty, the outermore ones are fully filled
            starting_row = nrows-len(rows)-1
            seats_on_starting_row = nseats-sum(rows)
            del rows

        case _:
            raise ValueError(f"Unrecognized strategy : {filling_strategy}")

    positions = {}
    for r in range(starting_row, nrows):
        if r == nrows-1: # if it's the last, outermost row
            # fit all the remaining seats
            nseats_this_row = nseats-len(positions)
        elif filling_strategy == FillingStrategy.OUTER_PRIORITY:
            if r == starting_row:
                nseats_this_row = seats_on_starting_row
            else:
                nseats_this_row = maxed_rows[r]
        else:
            # fullness of the diagram times the maximum number of seats in the row
            nseats_this_row = round(filling_ratio * maxed_rows[r])
            # actually more precise rounding : avoid rounding errors to accumulate too much
            # nseats_this_row = round((nseats-len(positions)) * maxed_rows[r]/sum(maxed_rows[r:]))

        # row radius : the radius of the circle crossing the center of each seat in the row
        row_arc_radius = .5 + 2*r*row_thicc

        if nseats_this_row == 1:
            positions[1., row_arc_radius] = math.pi/2
        else:
            # the angle necessary in this row to put the first (and last) seats fully in the canvas
            angle_margin = math.asin(row_thicc/row_arc_radius)
            # add the margin to make up the side angle
            angle_margin += span_angle_margin
            # alternatively, allow the centers of the seats by the side to reach the angle's boundary
            # angle_margin = max(angle_margin, span_angle_margin)

            # the angle separating two seats of that row
            angle_increment = (math.pi-2*angle_margin) / (nseats_this_row-1)
            # a fraction of the remaining space,
            # keeping in mind that the same elevation on start and end limits that remaining place to less than 2pi

            for s in range(nseats_this_row):
                angle = angle_margin + s*angle_increment
                # an oriented angle, so it goes trig positive (counterclockwise)
                positions[row_arc_radius*math.cos(angle)+1, row_arc_radius*math.sin(angle)] = angle

    return positions