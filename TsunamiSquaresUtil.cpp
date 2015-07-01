// Copyright (c) 2015 Kasey W. Schultz, Eric M. Heien
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "TsunamiSquaresUtil.h"

std::ostream &tsunamisquares::operator<<(std::ostream &os, const LatLonDepth &pt) {
    os << "[" << pt.lat() << "," << pt.lon();

    if (pt.altitude() != 0) os << "," << pt.altitude();

    os << "]";
    return os;
}

std::ostream &tsunamisquares::operator<<(std::ostream &os, const Vec<2> &v) {
    for (int i=0; i<2; ++i) os << (i == 0 ? "(" : ",") << v[i] << (i == 1 ? ")" : "");

    return os;
}

std::ostream &tsunamisquares::operator<<(std::ostream &os, const Vec<3> &v) {
    for (int i=0; i<3; ++i) os << (i == 0 ? "(" : ",") << v[i] << (i == 2 ? ")" : "");

    return os;
}

std::ostream &tsunamisquares::operator<<(std::ostream &os, const Tensor<3,3> &tensor) {
    for (int i=0; i<3; ++i) os << (i == 0 ? "[" : ",") << tensor[i] << (i == 2 ? "]" : "");

    return os;
}

std::ostream &tsunamisquares::operator<<(std::ostream &os, const TensorRow<3> &tr) {
    for (int i=0; i<3; ++i) os << (i == 0 ? "[" : ",") << tr[i] << (i == 2 ? "]" : "");

    return os;
}

/*!
 Converts an (x, y) value (given in km) to lattitude/longitude (in degrees decimal) relative to a set origin.
 This just inverts the conversion done in convert2xy
 */
#ifdef GEOGRAPHICLIB_FOUND
tsunamisquares::LatLonDepth tsunamisquares::Conversion::convert2LatLon(const Vec<3> &in_pt) const {
    double  a = GeographicLib::Constants::WGS84_a(),    // major radius
            f = GeographicLib::Constants::WGS84_f();    // flattening
    const   GeographicLib::Geodesic geod(a, f);
    double  s12, azi1, a12;
    double  new_lat, new_lon, new_altitude;

    s12 = sqrt( in_pt[0] * in_pt[0] + in_pt[1] * in_pt[1] );
    azi1 = rad2deg(atan2(in_pt[0],in_pt[1]));
    a12 = geod.Direct(_base.lat(),_base.lon(), azi1, s12, new_lat, new_lon);
    new_altitude = in_pt[2];

    return LatLonDepth(new_lat, new_lon, new_altitude);
}
#else
/*!
 Converts XYZ coordinate to latitude/longitude off _base coordinate using
 Vincenty inverse formula for ellipsoids.

 Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2011
 http://www.movable-type.co.uk/scripts/latlong.html

 from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the
 Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975
 http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
 */
tsunamisquares::LatLonDepth tsunamisquares::Conversion::convert2LatLon(const Vec<3> &in_pt) const {
    double a = EARTH_EQUATORIAL_RADIUS,  f = WGS_84_FLATTENING, b = (1-f)*a;
    double s = sqrt( in_pt[0] * in_pt[0] + in_pt[1] * in_pt[1] );
    double alpha1 = atan2(in_pt[0], in_pt[1]);
    double sinAlpha1 = sin(alpha1);
    double cosAlpha1 = cos(alpha1);

    double tanU1 = (1-f) * tan(deg2rad(_base.lat()));
    double cosU1 = 1.0 / sqrt((1 + tanU1*tanU1)), sinU1 = tanU1*cosU1;
    double sigma1 = atan2(tanU1, cosAlpha1);
    double sinAlpha = cosU1 * sinAlpha1;
    double cosSqAlpha = 1 - sinAlpha*sinAlpha;
    double uSq = cosSqAlpha * (a*a - b*b) / (b*b);
    double A = 1 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)));
    double B = uSq/1024.0 * (256.0+uSq*(-128.0+uSq*(74.0-47.0*uSq)));

    double sigma = s / (b*A), sigmaP = 2*M_PI;
    double cos2SigmaM, sinSigma, cosSigma, deltaSigma;

    while (fabs(sigma-sigmaP) > 1e-14) {
        cos2SigmaM = cos(2*sigma1 + sigma);
        sinSigma = sin(sigma);
        cosSigma = cos(sigma);
        deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
        sigmaP = sigma;
        sigma = s / (b*A) + deltaSigma;
    }

    double tmp = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
    double lat2 = atan2(sinU1*cosSigma + cosU1*sinSigma*cosAlpha1,
                        (1-f)*sqrt(sinAlpha*sinAlpha + tmp*tmp));
    double lambda = atan2(sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1);
    double C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
    double L = lambda - (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
    double lon2 = fmod((deg2rad(_base.lon())+L+3*M_PI),(2*M_PI)) - M_PI;  // normalise to -180...+180

    return LatLonDepth(rad2deg(lat2), rad2deg(lon2), in_pt[2]);
}
#endif
/*!
 Converts a latitude/longitude point to (x, y) in meters from a set origin.
 */
#ifdef GEOGRAPHICLIB_FOUND
tsunamisquares::Vec<3> tsunamisquares::Conversion::convert2xyz(const LatLonDepth &in_pt) const {
    double  a = GeographicLib::Constants::WGS84_a(),  // major radius
            f = GeographicLib::Constants::WGS84_f();  // flattening
    const   GeographicLib::Geodesic geod(a, f);
    double  new_vals[3];
    double  s12, azi1, azi2, a12;

    // s12 is the distance from (latitude0, longitude0) to (in_pt.lat,in_pt.lon)
    // azi1 is the direction (angle in degrees clockwise from north) the point (in_pt.lat,in_pt.lon) is in
    a12 = geod.Inverse(_base.lat(),_base.lon(),in_pt.lat(),in_pt.lon(),s12,azi1,azi2);

    new_vals[0] = s12 * sin(deg2rad(azi1));
    new_vals[1] = s12 * cos(deg2rad(azi1));
    new_vals[2] = in_pt.altitude();

    return Vec<3>(new_vals);
}
#else
tsunamisquares::Vec<3> tsunamisquares::Conversion::convert2xyz(const LatLonDepth &in_pt) const {
    double      new_vals[3];
    double      dist, start_azimuth, end_azimuth;

    dist_vincenty(dist, start_azimuth, end_azimuth, LatLonDepth(_base.lat(), _base.lon()), LatLonDepth(in_pt.lat(), in_pt.lon()));

    new_vals[0] = dist*sin(start_azimuth);
    new_vals[1] = dist*cos(start_azimuth);
    new_vals[2] = in_pt.altitude();

    return Vec<3>(new_vals);
}
#endif

tsunamisquares::VectorList tsunamisquares::Conversion::convertArray2xyz(const FloatList &lats, const FloatList &lons) const {
    tsunamisquares::VectorList conversions;

    for (FloatList::size_type lat_id = 0; lat_id != lats.size(); lat_id++) {
        for (FloatList::size_type lon_id = 0; lon_id != lons.size(); lon_id++) {
            conversions.push_back(convert2xyz(LatLonDepth(lats[lat_id],lons[lon_id],0.0)));
        }
    }

    return conversions;
};


/*!
 Calculates geodetic distance between two points specified by latitude/longitude using
 Vincenty inverse formula for ellipsoids.

 Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2011
 http://www.movable-type.co.uk/scripts/latlong.html

 from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the
       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975
       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
 */
void tsunamisquares::Conversion::dist_vincenty(double &distance, double &start_azimuth, double &end_azimuth, const LatLonDepth &p1, const LatLonDepth &p2) const {
    double      p1lon = deg2rad(p1.lon()), p1lat = deg2rad(p1.lat());
    double      p2lon = deg2rad(p2.lon()), p2lat = deg2rad(p2.lat());
    double      a = EARTH_EQUATORIAL_RADIUS;
    double      f = WGS_84_FLATTENING;
    double      b = (1-f)*a;
    double      L = p2lon-p1lon;
    double      U1 = atan((1-f)*tan(p1lat));
    double      U2 = atan((1-f)*tan(p2lat));
    double      sinU1 = sin(U1), cosU1 = cos(U1);
    double      sinU2 = sin(U2), cosU2 = cos(U2);
    double      lambda = L, lambdaP;
    int         max_iter = 100;
    double      sinLambda, cosLambda, sinSigma, cosSigma;
    double      sigma, sinAlpha, cosSqAlpha, cos2SigmaM;
    double      C, uSq, A, B, deltaSigma;

    do {
        sinLambda = sin(lambda);
        cosLambda = cos(lambda);
        sinSigma = sqrt(pow((cosU2*sinLambda),2) + pow((cosU1*sinU2-sinU1*cosU2*cosLambda),2));

        if (sinSigma == 0) {        // coincident points
            distance = start_azimuth = end_azimuth = 0;
            return;
        }

        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        sinAlpha = cosU1*cosU2*sinLambda/sinSigma;
        cosSqAlpha = 1.0 - sinAlpha*sinAlpha;
        cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;

        if (cosSqAlpha == 0) cos2SigmaM = 0;    // equatorial line: cosSqAlpha = 0

        C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
    } while (fabs(lambda-lambdaP) > 1e-14 && --max_iter>0);

    if (max_iter == 0) {
        distance = start_azimuth = end_azimuth = nan("");
        return;
    };

    uSq = cosSqAlpha * (a*a - b*b)/(b*b);

    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));

    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

    deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)- B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));

    distance = b*A*(sigma-deltaSigma);

    start_azimuth = atan2(cosU2*sinLambda, cosU1*sinU2-sinU1*cosU2*cosLambda);

    end_azimuth = atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
}