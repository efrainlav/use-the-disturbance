from __future__ import division, absolute_import
import numpy as np

def geodetic2cartesian(height, latitude, longitude,
                       major_semiaxis, minor_semiaxis):
    '''
    Transform geocentric geodetic coordinates into geocentric
    Cartesian coordinates.

    input

    height: numpy array 1D - vector containing the geometric height (in meters).
    latitude: numpy array 1D - vector containing the latitude (in degrees).
    longitude: numpy array 1D - vector containing the longitude (in degrees).
    major_semiaxis: float - major semiaxis of the reference
                    ellipsoid (in meters).
    minor_semiaxis: float - minor semiaxis of the reference
                    ellipsoid (in meters).

    output

    x: numpy array 1D - vector containing the x component of the Cartesian
       coordinates (in meters).
    y: numpy array 1D - vector containing the y component of the Cartesian
       coordinates (in meters).
    z: numpy array 1D - vector containing the z component of the Cartesian
       coordinates (in meters).
    '''

    h = np.asarray(height)
    lat = np.asarray(latitude)
    lon = np.asarray(longitude)

    assert (h.size == lat.size == lon.size), 'height, latitude \
and longitude must have the same number of elements'
    assert (major_semiaxis > minor_semiaxis), 'major_semiaxis must be greater \
than the minor_semiaxis'


    #Prime vertical radius of curvature
    N = prime_vertical_curv(major_semiaxis, minor_semiaxis, lat)

    # convert degrees to radians
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    aux = N + height

    # squared first eccentricity
    e2 = (major_semiaxis**2. - minor_semiaxis**2.)/(major_semiaxis**2.)

    clat = np.cos(lat)
    slat = np.sin(lat)
    clon = np.cos(lon)
    slon = np.sin(lon)

    x = aux*clat*clon
    y = aux*clat*slon
    z = (N*(1 - e2) + height)*slat

    return x, y, z

def cartesian2geodetic(X, Y, Z, major_semiaxis, minor_semiaxis, itmax = 5):
    '''
    Convert geocentric Cartesian coordinates into geocentric geodetic
    coordinates by using the Hirvonen-Moritz algorithm.

    input

    X: numpy array 1D or float - vector containing the x component of the
       Cartesian coordinates (in meters).
    Y: numpy array 1D or float - vector containing the y component of the
       Cartesian coordinates (in meters).
    Z: numpy array 1D or float - vector containing the z component of the
       Cartesian coordinates (in meters).
    major_semiaxis: float - major semiaxis of the reference
                    ellipsoid (in meters).
    minor_semiaxis: float - minor semiaxis of the reference
                    ellipsoid (in meters).
    itmax: int - maximum number of iterations in the Hirvonen-Moritz algorithm.

    output

    height: numpy array 1D - vector containing the geometric height (in meters).
    latitude: numpy array 1D - vector containing the latitude (in degrees).
    longitude: numpy array 1D - vector containing the longitude (in degrees).

    '''

    x = np.asarray(X)
    y = np.asarray(Y)
    z = np.asarray(Z)

    assert (x.size == y.size == z.size), \
'x, y and z must have the same number of elements'
    assert (major_semiaxis > minor_semiaxis), 'major_semiaxis must be greater \
than the minor_semiaxis'

    # horizontal distance
    p = np.sqrt(x**2. + y**2.)

    # null and non-null horizontal distances
    p_non_null = (p >= 1e-8)
    p_null = np.logical_not(p_non_null)

    lon = np.zeros_like(x)
    lat = np.zeros_like(x)
    height = np.zeros_like(x)

    # squared first eccentricity
    e2 = (major_semiaxis**2. - minor_semiaxis**2.)/(major_semiaxis**2.)

    aux1 = z[p_non_null]/p[p_non_null]
    aux2 = 1.- e2

    # define the coordinates for null horizontal distances
    lon[p_null] = 0.
    height[p_null] = np.abs(z[p_null]) - minor_semiaxis
    lat[p_null] = np.sign(z[p_null])*np.pi*0.5

    # first iteration
    lat[p_non_null] = np.arctan(aux1/aux2)
    sinlat = np.sin(lat[p_non_null])
    N = major_semiaxis/np.sqrt(1 - e2*sinlat*sinlat)
    height[p_non_null] = p[p_non_null]/np.cos(lat[p_non_null]) - N

    for i in range(itmax):
        aux3 = e2*N/(N + height[p_non_null])
        lat[p_non_null] = np.arctan(aux1/(1.-aux3))
        sinlat = np.sin(lat[p_non_null])
        N = major_semiaxis/np.sqrt(1 - e2*sinlat*sinlat)
        height[p_non_null] = p[p_non_null]/np.cos(lat[p_non_null]) - N

    lon[p_non_null] = np.arctan2(y[p_non_null], x[p_non_null])

    # convert latitude and longitude from radians to degrees
    latitude = np.rad2deg(lat)
    longitude = np.rad2deg(lon)

    return height, latitude, longitude


def cartesian2geodetic_approx(X, Y, Z, major_semiaxis, minor_semiaxis):
    '''
    Convert geocentric Cartesian coordinates into geocentric geodetic
    coordinates by using an approximated formula (Hofmann-Wellenhof and
    Moritz, 2005).

    reference

    Hofmann-Wellenhof, B. and H. Moritz, 2005, Physical Geodesy. Springer.

    input

    X: numpy array 1D or float - vector containing the x component of the
       Cartesian coordinates (in meters).
    Y: numpy array 1D or float - vector containing the y component of the
       Cartesian coordinates (in meters).
    Z: numpy array 1D or float - vector containing the z component of the
       Cartesian coordinates (in meters).
    major_semiaxis: float - major semiaxis of the reference
                    ellipsoid (in meters).
    minor_semiaxis: float - minor semiaxis of the reference
                    ellipsoid (in meters).
    itmax: int - maximum number of iterations in the Hirvonen-Moritz algorithm.

    output

    height: numpy array 1D - vector containing the geometric height (in meters).
    latitude: numpy array 1D - vector containing the latitude (in degrees).
    longitude: numpy array 1D - vector containing the longitude (in degrees).

    '''

    x = np.asarray(X)
    y = np.asarray(Y)
    z = np.asarray(Z)

    assert (x.size == y.size == z.size), \
'x, y and z must have the same number of elements'
    assert (major_semiaxis > minor_semiaxis), 'major_semiaxis must be greater \
than the minor_semiaxis'

    # horizontal distance
    p = np.sqrt(x**2. + y**2.)

    # null and non-null horizontal distances
    p_non_null = (p >= 1e-8)
    p_null = np.logical_not(p_non_null)

    lon = np.zeros_like(x)
    lat = np.zeros_like(x)
    height = np.zeros_like(x)

    # define the coordinates for null horizontal distances
    lon[p_null] = 0.
    height[p_null] = np.abs(z[p_null]) - minor_semiaxis
    lat[p_null] = np.sign(z[p_null])*np.pi*0.5

    # squared first eccentricity
    e2 = (major_semiaxis**2. - minor_semiaxis**2.)/(major_semiaxis**2.)

    # squared second eccentricity
    elinha2 = (major_semiaxis**2. - minor_semiaxis**2.)/(minor_semiaxis**2.)

    # auxiliary variable
    theta = np.arctan(z[p_non_null]*major_semiaxis/\
                      (p[p_non_null]*minor_semiaxis))
    sintheta = np.sin(theta)
    costheta = np.cos(theta)

    aux1 = z[p_non_null] + elinha2*minor_semiaxis*sintheta*sintheta*sintheta
    aux2 = p[p_non_null] - e2*major_semiaxis*costheta*costheta*costheta

    #lat[p_non_null] = np.arctan(aux1/aux2)
    lat[p_non_null] = np.arctan2(aux1, aux2)
    #lon[p_non_null] = np.arctan(y[p_non_null]/x[p_non_null])
    lon[p_non_null] = np.arctan2(y[p_non_null], x[p_non_null])

    sinlat = np.sin(lat[p_non_null])
    N = major_semiaxis/np.sqrt(1 - e2*sinlat*sinlat)

    height[p_non_null] = p[p_non_null]/np.cos(lat[p_non_null]) - N

    # convert latitude and longitude from radians to degrees
    latitude = np.rad2deg(lat)
    longitude = np.rad2deg(lon)

    return height, latitude, longitude


def molodensky_completa(h1, phi1, lamb1, a1, f1, a2, f2, dx, dy, dz):
    '''
    Transforma coordenadas geodesicas h1 (altitude), phi1 (latitude) e
    lamb1 (longitude) - referidas a um elipsoide 1 - em coordenadas geodesicas
    h2 (altitude), phi2 (latitude) e lamb2 (longitude) - referidas a um
    elipsoide 2.

    input

    h1: numpy array 1D - vetor de altitudes no elipsoide 1 (em metros).
    phi1: numpy array 1D - vetor de latitudes no elipsoide 1 (em radianos).
    lamb1: numpy array 1D - vetor de longitudes no elipsoide 1 (em radianos).
    a1: float - semieixo maior do elipsoide 1.
    f1: float - achatamento do elipsoide 1.
    a2: float - semieixo maior do elipsoide 2.
    f2: float - achatamento do elipsoide 2.
    dx: float - translacao, ao longo do eixo x, da origem do elipsoide 2
        em relacao ao elipsoide 1.
    dy: float - translacao, ao longo do eixo y, da origem do elipsoide 2
        em relacao ao elipsoide 1.
    dz: float - translacao, ao longo do eixo z, da origem do elipsoide 2
        em relacao ao elipsoide 1.

    output

    h2: numpy array 1D - vetor de altitudes no elipsoide 2 (em metros).
    phi2: numpy array 1D - vetor de latitudes no elipsoide 2 (em radianos).
    lamb2: numpy array 1D - vetor de longitudes no elipsoide 2 (em radianos).
    '''

    assert (h1.size == phi1.size == lamb1.size), \
        'h1, phi1 e lamb1 devem ter o mesmo numero de elementos'
    assert (f1 < a1), 'f1 deve ser menor que a1'
    assert (f2 < a2), 'f2 deve ser menor que a2'

    senophi = np.sin(phi1)
    cossenophi = np.cos(phi1)

    senolamb = np.sin(lamb1)
    cossenolamb = np.cos(lamb1)

    da = a2 - a1
    df = f2 - f1
    a = a1 + 0.5*da
    f = f1 + 0.5*df
    b = (1. - f)*a
    e2 = (a**2. - b**2.)/a**2.

    W = np.sqrt(1. - e2*(senophi**2.))
    W[W < 1.0e-8] == 1.0e-8
    M = a*(1. - e2)/(W**3.)
    N = a/W

    m11 = -senophi*cossenolamb
    m12 = -senophi*senolamb
    m13 = cossenophi
    m14 = e2*senophi*cossenophi/W
    m15 = senophi*cossenophi*(M*a/b + N*b/a)
    m21 = -senolamb
    m22 = cossenolamb
    m31 = cossenophi*cossenolamb
    m32 = cossenophi*senolamb
    m33 = senophi
    m34 = -W
    m35 = a*(1. - f)*(senophi**2.)/W

    dphi = (m11*dx + m12*dy + m13*dz + m14*da + m15*df)/(M + h1)
    dlamb = (m21*dx + m22*dy)/((N + h1)*cossenophi)
    dh = m31*dx + m32*dy + m33*dz + m34*da + m35*df

    h2 = h1 + dh
    phi2 = phi1 + dphi
    lamb2 = lamb1 + dlamb

    return h2, phi2, lamb2

def prime_vertical_curv(a, b, latitude):
    '''
    Compute the prime vertical radius of curvature.

    input

    a: float - major semiaxis (in meters) of the reference ellipsoid.
    b: float - minor semiaxis (in meters) of the reference ellipsoid.
    latitude: numpy array 1D - latitude (in degrees) of the computation points.

    output

    N: numpy array 1D - prime vertical radius of curvature
       computed at each latitude.
    '''

    assert a > b, 'major semiaxis must be greater than minro semiaxis'
    assert a > 0, 'major semiaxis must be nonnull'
    assert b > 0, 'minor semiaxis must be nonnull'

    # squared first eccentricity
    e2 = (a*a - b*b)/(a*a)

    # squared sine of latitude
    sin2lat = np.sin(np.deg2rad(latitude))
    sin2lat *= sin2lat

    N = a/np.sqrt(1. - e2*sin2lat)

    return N

def meridian_curv(a, b, latitude):
    '''
    Compute the prime vertical radius of curvature.

    input

    a: float - major semiaxis (in meters) of the reference ellipsoid.
    b: float - minor semiaxis (in meters) of the reference ellipsoid.
    latitude: numpy array 1D - latitude (in degrees) of the computation points.

    output

    M: numpy array 1D - meridian radius of curvature
       computed at each latitude.
    '''

    assert a > b, 'major semiaxis must be greater than minro semiaxis'
    assert a > 0, 'major semiaxis must be nonnull'
    assert b > 0, 'minor semiaxis must be nonnull'

    # squared first eccentricity
    e2 = (a*a - b*b)/(a*a)

    # squared sine of latitude
    sin2lat = np.sin(np.deg2rad(latitude))
    sin2lat *= sin2lat

    # auxiliary variable
    aux = np.sqrt(1. - e2*sin2lat)

    M = (a*(1 - e2))/(aux*aux*aux)

    return M


def rotation_matrix(latitude, longitude):
    '''
    Compute the elements of a rotation matrix whose columns are
    unit vectors pointing to the direction of increasing orthometric height,
    the direction of increasing geodetic latitude and the direction of
    increasing longitude. At a given point with geodetic coordinates
    (h, lat, lon), the rotation matrix can be written as follows:

        | cos(lat)*cos(lon)   -sin(lat)*cos(lon)   -sin(lon) |
    R = | cos(lat)*sin(lon)   -sin(lat)*sin(lon)    cos(lon) |
        | sin(lat)             cos(lat)             0        |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    R: list of numpy arrays 1D - list of vectors containing the elements
       11, 12, 13, 21, 22, 23, 31 and 32 of the rotation matrix evaluated
       at the computation points.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    R11 = coslat*coslon
    R12 = -sinlat*coslon
    R13 = -sinlon
    R21 = coslat*sinlon
    R22 = -sinlat*sinlon
    R23 = coslon
    R31 = sinlat
    R32 = coslat

    R = [R11, R12, R13, R21, R22, R23, R31, R32]

    return R


def rotation_matrix(latitude, longitude):
    '''
    Compute the elements of a rotation matrix whose columns are
    unit vectors pointing to the direction of increasing orthometric height,
    the direction of increasing geodetic latitude and the direction of
    increasing longitude. At a given point with geodetic coordinates
    (h, lat, lon), the rotation matrix can be written as follows:

        | cos(lat)*cos(lon)   -sin(lat)*cos(lon)   -sin(lon) |
    R = | cos(lat)*sin(lon)   -sin(lat)*sin(lon)    cos(lon) |
        | sin(lat)             cos(lat)             0        |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    R: list of numpy arrays 1D - list of vectors containing the elements
       11, 12, 13, 21, 22, 23, 31 and 32 of the rotation matrix evaluated
       at the computation points.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    R11 = coslat*coslon
    R12 = -sinlat*coslon
    R13 = -sinlon
    R21 = coslat*sinlon
    R22 = -sinlat*sinlon
    R23 = coslon
    R31 = sinlat
    R32 = coslat

    R = [R11, R12, R13, R21, R22, R23, R31, R32]

    return R
