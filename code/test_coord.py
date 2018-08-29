from __future__ import division, absolute_import
import numpy as np
from code.cartography import coord
from numpy.testing import assert_almost_equal
from pytest import raises


def test_prime_radius_bad_semiaxes():
    'major semiaxis lower than minor semiaxis'
    latitude = np.ones(100)
    raises(AssertionError, coord.prime_vertical_curv, 10, 23, latitude)


def test_prime_radius_negative_semiaxes():
    'null semiaxis'
    latitude = np.ones(100)
    raises(AssertionError, coord.prime_vertical_curv, -10, -23, latitude)


def test_meridian_radius_bad_semiaxes():
    'major semiaxis lower than minor semiaxis'
    latitude = np.ones(100)
    raises(AssertionError, coord.meridian_curv, 10, 23, latitude)


def test_meridian_radius_negative_semiaxes():
    'null semiaxis'
    latitude = np.ones(100)
    raises(AssertionError, coord.meridian_curv, -10, -23, latitude)


def test_relationship_curvatures():
    'verify relationship between the curvatures'
    latitude = np.linspace(-90, 90, 181)
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)
    N = coord.prime_vertical_curv(a, b, latitude)
    M = coord.meridian_curv(a, b, latitude)
    e2 = (a*a - b*b)/(a*a)
    sin2lat = np.sin(np.deg2rad(latitude))
    sin2lat *= sin2lat
    M_relationship = ((1 - e2)/(1 - e2*sin2lat))*N
    assert_almost_equal(M, M_relationship, decimal=8)


def test_prime_radius_known_input():
    'verify results obtained for known input'
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)
    e2 = (a*a - b*b)/(a*a)
    N_true_0 = a
    N_true_90 = a/np.sqrt(1 - e2)
    N_calc_0 = coord.prime_vertical_curv(a, b, 0)
    N_calc_90 = coord.prime_vertical_curv(a, b, 90)
    assert_almost_equal(N_true_0, N_calc_0, decimal=15)
    assert_almost_equal(N_true_90, N_calc_90, decimal=15)


def test_prime_radius_known_input_2():
    'verify results obtained for known input'
    # true geodetic coordinates
    lat_true = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # squared first eccentricity
    e2 = (a**2. - b**2.)/(a**2.)

    # true prime vertical radius of curvature
    N_true = a/np.sqrt(1 - e2*sinlat*sinlat)

    # computed prime vertical radius of curvature
    N = coord.prime_vertical_curv(a, b, lat_true)

    assert_almost_equal(N_true, N, decimal=15)


def test_meridian_radius_known_input():
    'verify results obtained for known input'
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)
    e2 = (a*a - b*b)/(a*a)
    M_true_0 = a*(1 - e2)
    M_true_90 = a/np.sqrt(1 - e2)
    M_calc_0 = coord.meridian_curv(a, b, 0)
    M_calc_90 = coord.meridian_curv(a, b, 90)
    assert_almost_equal(M_true_0, M_calc_0, decimal=15)
    assert_almost_equal(M_true_90, M_calc_90, decimal=8)


def test_meridian_radius_known_input_2():
    'verify results obtained for known input'
    # true geodetic coordinates
    lat_true = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # squared first eccentricity
    e2 = (a**2. - b**2.)/(a**2.)

    # true meridian radius of curvature
    M_true = (a*(1 - e2))/np.power(1 - e2*sinlat*sinlat, 1.5)

    # computed meridian radius of curvature
    M = coord.meridian_curv(a, b, lat_true)

    assert_almost_equal(M_true, M, decimal=8)


def test_rotation_matrix_known_values():
    'verify results obtained for known input'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([0, 0, 0, 90, 90, 90, 180, 180])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])
    coslat = np.array([1, (np.sqrt(6) + np.sqrt(2))/4,
                      (np.sqrt(2 + np.sqrt(2)))/2, np.sqrt(3)/2,
                      np.sqrt(2)/2, 0.5, (np.sqrt(6) - np.sqrt(2))/4, 0])
    sinlon = np.array([0, 0, 0, 1, 1, 1, 0, 0])
    coslon = np.array([1, 1, 1, 0, 0, 0, -1, -1])

    R11_true = coslat*coslon
    R12_true = -sinlat*coslon
    R13_true = -sinlon
    R21_true = coslat*sinlon
    R22_true = -sinlat*sinlon
    R23_true = coslon
    R31_true = sinlat
    R32_true = coslat

    R = coord.rotation_matrix(latitude, longitude)

    assert_almost_equal(R11_true, R[0], decimal=15)
    assert_almost_equal(R12_true, R[1], decimal=15)
    assert_almost_equal(R13_true, R[2], decimal=15)
    assert_almost_equal(R21_true, R[3], decimal=15)
    assert_almost_equal(R22_true, R[4], decimal=15)
    assert_almost_equal(R23_true, R[5], decimal=15)
    assert_almost_equal(R31_true, R[6], decimal=15)
    assert_almost_equal(R32_true, R[7], decimal=15)


def test_rotation_matrix_known_values_floats():
    'verify results obtained for known input'
    latitude = -15
    longitude = 0

    sinlat = -(np.sqrt(6) - np.sqrt(2))/4
    coslat = (np.sqrt(6) + np.sqrt(2))/4
    sinlon = 0
    coslon = 1

    R11_true = coslat*coslon
    R12_true = -sinlat*coslon
    R13_true = -sinlon
    R21_true = coslat*sinlon
    R22_true = -sinlat*sinlon
    R23_true = coslon
    R31_true = sinlat
    R32_true = coslat

    R = coord.rotation_matrix(latitude, longitude)

    assert_almost_equal(R11_true, R[0], decimal=15)
    assert_almost_equal(R12_true, R[1], decimal=15)
    assert_almost_equal(R13_true, R[2], decimal=15)
    assert_almost_equal(R21_true, R[3], decimal=15)
    assert_almost_equal(R22_true, R[4], decimal=15)
    assert_almost_equal(R23_true, R[5], decimal=15)
    assert_almost_equal(R31_true, R[6], decimal=15)
    assert_almost_equal(R32_true, R[7], decimal=15)


def test_rotation_matrix_orthogonality():
    'Rotation matrix must be mutually orthogonal'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([-17, 0, 30, 9, 90, 23, 180, 1])
    R = coord.rotation_matrix(latitude, longitude)
    for R11i, R12i, R13i, R21i, R22i, R23i, R31i, R32i in \
        zip(R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7]):
        Ri = np.array([[R11i, R12i, R13i],
                       [R21i, R22i, R23i],
                       [R31i, R32i, 0]])
        assert_almost_equal(np.dot(Ri.T, Ri), np.identity(3), decimal=15)
        assert_almost_equal(np.dot(Ri, Ri.T), np.identity(3), decimal=15)


def test_rotation_matrix_lines_bad_arguments():
    'latitude and longitude with different number of elements'
    latitude = np.ones(100)
    longitude = np.zeros(34)
    raises(AssertionError, coord.rotation_matrix, latitude, longitude)


def test_geodetic2cartesian_bad_coordinates():
    'heigth, latitude and longitude with different number of elements'

    major_semiaxis = 10
    minor_semiaxes = 8

    latitude = np.empty(100)
    longitude = np.empty(34)
    height = np.empty(34)
    raises(AssertionError, coord.geodetic2cartesian, height, latitude,
           longitude, major_semiaxis, minor_semiaxes)

    latitude = np.empty(34)
    longitude = np.empty(100)
    height = np.empty(34)
    raises(AssertionError, coord.geodetic2cartesian, height, latitude,
           longitude, major_semiaxis, minor_semiaxes)

    latitude = np.empty(34)
    longitude = np.empty(34)
    height = np.empty(100)
    raises(AssertionError, coord.geodetic2cartesian, height, latitude,
           longitude, major_semiaxis, minor_semiaxes)


def test_geodetic2cartesian_bad_semiaxes():
    'major semiaxis lower than minor semiaxis'

    major_semiaxis = 17
    minor_semiaxes = 20

    latitude = np.empty(34)
    longitude = np.empty(34)
    height = np.empty(34)
    raises(AssertionError, coord.geodetic2cartesian, height, latitude,
           longitude, major_semiaxis, minor_semiaxes)


def test_geodetic2cartesian_known_input():
    'verify the computed coordinates obtained from specific input'

    # geometric height (meters), latitude (degrees) and longitude (degrees)
    h = np.array([0, 0, 0])
    lat = np.array([0, 90, 0])
    lon = np.array([0, 0, 90])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # true coordinates
    x_true = np.array([a, 0, 0])
    y_true = np.array([0, 0, a])
    z_true = np.array([0, b, 0])

    # computed coordinates
    x, y, z = coord.geodetic2cartesian(h, lat, lon, a, b)

    assert_almost_equal(x_true, x, decimal=9)
    assert_almost_equal(y_true, y, decimal=15)
    assert_almost_equal(z_true, z, decimal=15)


def test_geodetic2cartesian_known_input_2():
    'verify the computed coordinates obtained from specific input'

    # true geodetic coordinates
    h_true = np.array([0, -100, 200, -300, 400, -800, 1600, 7000])
    lat_true = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    lon_true = np.array([0, 0, 0, 90, 90, 90, 180, 180])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])
    coslat = np.array([1, (np.sqrt(6) + np.sqrt(2))/4,
                      (np.sqrt(2 + np.sqrt(2)))/2, np.sqrt(3)/2,
                      np.sqrt(2)/2, 0.5, (np.sqrt(6) - np.sqrt(2))/4, 0])
    sinlon = np.array([0, 0, 0, 1, 1, 1, 0, 0])
    coslon = np.array([1, 1, 1, 0, 0, 0, -1, -1])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # squared first eccentricity
    e2 = (a**2. - b**2.)/(a**2.)

    # prime vertical radius of curvature
    N = coord.prime_vertical_curv(a, b, lat_true)

    # true Cartesian coordinates
    x_true = (N + h_true)*coslat*coslon
    y_true = (N + h_true)*coslat*sinlon
    z_true = (N*(1 - e2) + h_true)*sinlat

    # computed Cartesian coordinates
    x, y, z = coord.geodetic2cartesian(h_true, lat_true, lon_true, a, b)

    assert_almost_equal(x_true, x, decimal=9)
    assert_almost_equal(y_true, y, decimal=9)
    assert_almost_equal(z_true, z, decimal=9)


def test_cartesian2geodetic_bad_coordinates():
    'x, y and z with different number of elements'

    major_semiaxis = 10
    minor_semiaxes = 8

    x = np.empty(100)
    y = np.empty(34)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic, x, y, z,
           major_semiaxis, minor_semiaxes)

    x = np.empty(34)
    y = np.empty(100)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic, x, y, z,
           major_semiaxis, minor_semiaxes)

    x = np.empty(34)
    y = np.empty(34)
    z = np.empty(100)
    raises(AssertionError, coord.cartesian2geodetic, x, y, z,
           major_semiaxis, minor_semiaxes)


def test_cartesian2geodetic_bad_semiaxes():
    'major semiaxis lower than minor semiaxis'

    major_semiaxis = 17
    minor_semiaxes = 20

    x = np.empty(34)
    y = np.empty(34)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic, x, y, z,
           major_semiaxis, minor_semiaxes)


def test_cartesian2geodetic_known_input():
    'verify the computed coordinates obtained from specific input'

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # Cartesian components x, y and z (in meters)
    x = np.array([a + 700, 0, 0])
    y = np.array([0, 0, a - 1200])
    z = np.array([0, b, 0])

    # true coordinates
    h_true = np.array([700, 0, -1200])
    lat_true = np.array([0, 90, 0])
    lon_true = np.array([0, 0, 90])

    # computed coordinates
    h, lat, lon = coord.cartesian2geodetic(x, y, z, a, b)

    assert_almost_equal(h_true, h, decimal=15)
    assert_almost_equal(lat_true, lat, decimal=15)
    assert_almost_equal(lon_true, lon, decimal=15)


def test_cartesian2geodetic_known_input_2():
    'verify the computed coordinates obtained from specific input'

    # true geodetic coordinates
    h_true = np.array([0, -100, 200, -300, 400, -800, 1600, 7000])
    lat_true = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    lon_true = np.array([0, 0, 0, 90, 90, 90, 180, 0])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])
    coslat = np.array([1, (np.sqrt(6) + np.sqrt(2))/4,
                      (np.sqrt(2 + np.sqrt(2)))/2, np.sqrt(3)/2,
                      np.sqrt(2)/2, 0.5, (np.sqrt(6) - np.sqrt(2))/4, 0])
    sinlon = np.array([0, 0, 0, 1, 1, 1, 0, 0])
    coslon = np.array([1, 1, 1, 0, 0, 0, -1, -1])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # squared first eccentricity
    e2 = (a**2. - b**2.)/(a**2.)

    # prime vertical radius of curvature
    N = coord.prime_vertical_curv(a, b, lat_true)

    # true Cartesian coordinates
    x_true = (N + h_true)*coslat*coslon
    y_true = (N + h_true)*coslat*sinlon
    z_true = (N*(1 - e2) + h_true)*sinlat

    # computed geodetic coordinates
    h, lat, lon = coord.cartesian2geodetic(x_true, y_true, z_true, a, b)

    assert_almost_equal(lat_true, lat, decimal=14)
    assert_almost_equal(lon_true, lon, decimal=15)
    assert_almost_equal(h_true, h, decimal=8)


def test_cartesian2geodetic_approx_bad_coordinates():
    'x, y and z with different number of elements'

    major_semiaxis = 10
    minor_semiaxes = 8

    x = np.empty(100)
    y = np.empty(34)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic_approx, x, y, z,
           major_semiaxis, minor_semiaxes)

    x = np.empty(34)
    y = np.empty(100)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic_approx, x, y, z,
           major_semiaxis, minor_semiaxes)

    x = np.empty(34)
    y = np.empty(34)
    z = np.empty(100)
    raises(AssertionError, coord.cartesian2geodetic_approx, x, y, z,
           major_semiaxis, minor_semiaxes)


def test_cartesian2geodetic_approx_bad_semiaxes():
    'major semiaxis lower than minor semiaxis'

    major_semiaxis = 17
    minor_semiaxes = 20

    x = np.empty(34)
    y = np.empty(34)
    z = np.empty(34)
    raises(AssertionError, coord.cartesian2geodetic_approx, x, y, z,
           major_semiaxis, minor_semiaxes)


def test_cartesian2geodetic_approx_known_input():
    'verify the computed coordinates obtained from specific input'

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # Cartesian components x, y and z (in meters)
    x = np.array([a + 700, 0, 0])
    y = np.array([0, 0, a - 1200])
    z = np.array([0, b, 0])

    # true coordinates
    h_true = np.array([700, 0, -1200])
    lat_true = np.array([0, 90, 0])
    lon_true = np.array([0, 0, 90])

    # computed coordinates
    h, lat, lon = coord.cartesian2geodetic_approx(x, y, z, a, b)

    assert_almost_equal(h_true, h, decimal=15)
    assert_almost_equal(lat_true, lat, decimal=15)
    assert_almost_equal(lon_true, lon, decimal=15)


def test_cartesian2geodetic_approx_known_input_2():
    'verify the computed coordinates obtained from specific input'

    # true geodetic coordinates
    lat_true = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    lon_true = np.array([0, 0, 0, 90, 90, 90, 180, 0])
    #h_true = np.array([0, -100, 200, -300, 400, -800, 1600, 7000])
    h_true = np.zeros_like(lat_true)

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])
    coslat = np.array([1, (np.sqrt(6) + np.sqrt(2))/4,
                      (np.sqrt(2 + np.sqrt(2)))/2, np.sqrt(3)/2,
                      np.sqrt(2)/2, 0.5, (np.sqrt(6) - np.sqrt(2))/4, 0])
    sinlon = np.array([0, 0, 0, 1, 1, 1, 0, 0])
    coslon = np.array([1, 1, 1, 0, 0, 0, -1, -1])

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # squared first eccentricity
    e2 = (a**2. - b**2.)/(a**2.)

    # prime vertical radius of curvature
    N = coord.prime_vertical_curv(a, b, lat_true)

    # true Cartesian coordinates
    x_true = (N + h_true)*coslat*coslon
    y_true = (N + h_true)*coslat*sinlon
    z_true = (N*(1 - e2) + h_true)*sinlat

    # computed geodetic coordinates
    h, lat, lon = coord.cartesian2geodetic_approx(x_true, y_true, z_true, a, b)

    assert_almost_equal(lat_true, lat, decimal=13)
    assert_almost_equal(lon_true, lon, decimal=15)
    assert_almost_equal(h_true, h, decimal=8)


def test_cartesian2geodetic_versus_geodetic2cartesian():
    'compare the results produced by cartesian2geodetic and geodetic2cartesian'

    # true geodetic coordinates
    h_true = np.linspace(-1000, 1000, 10)
    lat_true = np.linspace(-90, 89, 10)
    lon_true = np.linspace(0, 180, 10)

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # computed Cartesian coordinates
    x, y, z = coord.geodetic2cartesian(h_true, lat_true, lon_true, a, b)

    # computed geodetic coordinates
    h, lat, lon = coord.cartesian2geodetic(x, y, z, a, b)

    assert_almost_equal(h_true, h, decimal=8)
    assert_almost_equal(lat_true, lat, decimal=14)
    assert_almost_equal(lon_true, lon, decimal=14)


def test_cartesian2geodetic_versus_geodetic2cartesian_approx():
    'compare the results produced by cartesian2geodetic and \
geodetic2cartesian_approx'

    # true geodetic coordinates
    h_true = np.linspace(-1000, 1000, 10)
    lat_true = np.linspace(-90, 89, 10)
    lon_true = np.linspace(0, 180, 10)

    # major semiaxis, flattening, minor semiaxis and squared first eccentricity
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1-f)

    # computed Cartesian coordinates
    x, y, z = coord.geodetic2cartesian(h_true, lat_true, lon_true, a, b)

    # computed geodetic coordinates
    h, lat, lon = coord.cartesian2geodetic_approx(x, y, z, a, b)

    assert_almost_equal(h_true, h, decimal=8)
    assert_almost_equal(lat_true, lat, decimal=13)
    assert_almost_equal(lon_true, lon, decimal=14)
