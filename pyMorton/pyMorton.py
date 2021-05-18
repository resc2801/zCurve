from typing import List, Tuple
from gmpy2 import xmpz
import multiprocessing as mp


def interlace(*data_point: int, dims: int = None, bits_per_dim: int = None) -> int:
    """
    Interlace a given multi-dimensional data point into its 1D Morton code point.

    :param data_point:
         A multi-dimensional data point of unspecified dimensionality to encode.
    :param dims: int, optional
        The dimensionality of the underlying data space.; will speed things up if given.
    :param bits_per_dim: int, optional
        The number of encoding bits per dimension; will speed things up if given.
    :return: int
        A 1D Morton code point
    """
    dims = len(data_point) if dims is None else dims
    bits_per_dim = max(1, *data_point).bit_length() if bits_per_dim is None else bits_per_dim
    #    bits_per_dim = max(1, max(data_point).bit_length()) if bits_per_dim is None else bits_per_dim
    total_bits = dims * bits_per_dim
    c = xmpz()

    for i, v in enumerate(data_point):
        c[i:total_bits:dims] = v

    return int(c)


def par_interlace(data_points: List[List[int]], dims: int = None, bits_per_dim: int = None) -> List[int]:
    """

    :param data_points:
    :param dims:
    :param bits_per_dim:
    :return:
    """
    from functools import partial

    _interlace = partial(interlace, dims=dims, bits_per_dim=bits_per_dim)

    with mp.Pool(mp.cpu_count()) as pool:
        code_points = pool.starmap(_interlace, data_points, chunksize=None)

    return code_points


def deinterlace(code_point: int, dims: int = 3) -> List[int]:
    """
    Deinterlace given 1D Morton code_point into a multi-dimensional data point.

    :param code_point: int
        A 1D Morton code point
    :param dims: int
        The dimensionality of the underlying data space.
    :return: List[int]
        A multi-dimensional data point.
    """
    total_bits = code_point.bit_length() + (dims - code_point.bit_length() % dims)
    c = xmpz(code_point)
    data_point = [None] * dims

    for i in range(0, dims):
        data_point[i] = int(c[i:total_bits:dims])

    return data_point


def prev_morton(code_point: int, rmin_code: int, rmax_code: int, dims: int = 3) -> int:
    """
    Return 1D Morton code point previous to given 1D Morton code_point within range [rmin_code, rmax_code].

    :param code_point: int
        A 1D Morton code point
    :param rmin_code: int
        The minimum range 1D Morton code point
    :param rmax_code: int
        The maximum range 1D Morton code point
    :param dims:
        The dimensionality of the underlying data space.
    :return: int
        The 1D Morton code point previous to code_point within range [rmin_code, rmax_code].
    """
    total_bits = max(1, code_point, rmin_code, rmax_code).bit_length()
    total_bits = total_bits + (dims - total_bits % dims)

    CODE, MIN, MAX, LITMAX = (xmpz(code_point), xmpz(rmin_code), xmpz(rmax_code), xmpz())

    # bitwise scanning the codes of code_point and range_min_code and range_max_code starting from MSB
    for i in reversed(range(0, total_bits)):

        # the three bits are examined according to LITMAX decision table
        CODE_BIT, MIN_BIT, MAX_BIT = (CODE[i], MIN[i], MAX[i])

        if not CODE_BIT and not MIN_BIT and not MAX_BIT:
            # No action; continue
            continue
        elif not CODE_BIT and not MIN_BIT and MAX_BIT:
            # MAX = LOAD("0111...", MAX)
            MAX[i % dims:i:dims] = 2 ** (i + 1) - 1
            MAX[i] = 0
            continue
        elif not CODE_BIT and MIN_BIT and MAX_BIT:
            # finish
            return int(LITMAX)
        elif CODE_BIT and not MIN_BIT and not MAX_BIT:
            # LITMAX = MAX; finish
            return int(MAX)
        elif CODE_BIT and not MIN_BIT and MAX_BIT:
            # LITMAX = LOAD("0111...", MAX)
            LITMAX = MAX.copy()
            LITMAX[i % dims:i:dims] = 2**(i+1) - 1
            LITMAX[i] = 0
            # MIN = LOAD("10000...", MIN)
            MIN[i % dims:i:dims] = 0
            MIN[i] = 1
            continue
        elif CODE_BIT and MIN_BIT and MAX_BIT:
            # No action; continue
            continue
        else:
            raise ValueError("This case not possible because MIN <= MAX")

    return int(LITMAX)


def next_morton(code_point: int, rmin_code: int, rmax_code: int, dims: int = 3) -> int:
    """
    Return 1D Morton code point next to given 1D Morton code_point within range [rmin_code, rmax_code].

    :param code_point: int
        A 1D Morton code point
    :param rmin_code: int
        The minimum range 1D Morton code point
    :param rmax_code: int
        The maximum range 1D Morton code point
    :param dims:
        The dimensionality of the underlying data space.
    :return: int
        The 1D Morton code point next to code_point within range [rmin_code, rmax_code].
    """
    total_bits = max(1, code_point, rmin_code, rmax_code).bit_length()
    total_bits = total_bits + (dims - total_bits % dims)

    CODE, MIN, MAX, BIGMIN = (xmpz(code_point), xmpz(rmin_code), xmpz(rmax_code), xmpz())

    # bitwise scanning the codes of code_point and range_min_code and range_max_code starting from MSB
    for i in reversed(range(0, total_bits)):

        # the three bits are examined according to BIGMIN decision table
        CODE_BIT, MIN_BIT, MAX_BIT = (CODE[i], MIN[i], MAX[i])

        if not CODE_BIT and not MIN_BIT and not MAX_BIT:
            # No action; continue
            continue
        elif not CODE_BIT and not MIN_BIT and MAX_BIT:
            # BIGMIN=LOAD("1000...". MIN)
            BIGMIN = MIN.copy()
            BIGMIN[i % dims:i:dims] = 0
            BIGMIN[i] = 1
            # MAX = LOAD("0111...", MAX)
            MAX[i % dims:i:dims] = 2**(i+1) - 1
            MAX[i] = 0
            continue
        elif not CODE_BIT and MIN_BIT and MAX_BIT:
            # BIGMIN = MIN; finish
            return int(MIN)
        elif CODE_BIT and not MIN_BIT and not MAX_BIT:
            # finish
            return int(BIGMIN)
        elif CODE_BIT and not MIN_BIT and MAX_BIT:
            # MIN = LOAD("10000...", MIN)
            MIN[i % dims:i:dims] = 0
            MIN[i] = 1
            continue
        elif CODE_BIT and MIN_BIT and MAX_BIT:
            # No action; continue
            continue
        else:
            raise ValueError("This case not possible because MIN <= MAX")

    return int(BIGMIN)


def in_range(code_point: int, rmin_code: int, rmax_code: int, dims: int = 3) -> bool:
    if code_point < rmin_code or rmax_code < code_point:
        return False
    elif code_point == rmin_code or code_point == rmax_code:
        return True
    else:
        return code_point == prev_morton(next_morton(code_point, rmin_code, rmax_code, dims), rmin_code, rmax_code, dims)

