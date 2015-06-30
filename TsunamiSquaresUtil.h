// Copyright (c) 2012-2014 Eric M. Heien, Kasey W. Schultz
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

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <limits>

#define assertThrow(COND, ERR_MSG) assert(COND);

#ifdef GEOGRAPHICLIB_FOUND
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#endif

#define EARTH_MEAN_RADIUS       6371000.0           // in m
#define EARTH_EQUATORIAL_RADIUS 6378137.0           // in m
#define EARTH_POLAR_RADIUS      6356752.314245      // in m
#define WGS_84_FLATTENING       1/298.257223563

#define CROSS_TOLERANCE 0.000001

namespace tsunamisquares {
    // Provides general setup information regarding the library.
    std::string SetupInfo(void);

    /*!
     Represents a vector with associated operations.
     */
    template <unsigned int dim>
    class Vec {
        private:
            double _x[dim];

            void rotation_matrix(double rot_matrix[dim *dim], const double &theta) throw(std::domain_error) {
                if (dim != 3) {
                    throw std::domain_error("Vec::rotation_matrix");
                } else {
                    //[r11,r12,r13,r21,r22,r23,r31,r32,r33]
                    Vec<dim>        rot_axis;
                    double          w[dim];
                    unsigned int    i;

                    rot_axis = this->unit_vector();

                    for (i=0; i<dim; ++i) w[i] = rot_axis[i];

                    rot_matrix[0] = cos(theta) + pow(w[0],2.0) * (1.0 - cos(theta));
                    rot_matrix[1] = w[0] * w[1] * (1.0 - cos(theta)) + w[2] * sin(theta);
                    rot_matrix[2] = w[0] * w[2] * (1.0 - cos(theta)) - w[1] * sin(theta);
                    rot_matrix[3] = w[1] * w[0] * (1.0 - cos(theta)) - w[2] * sin(theta);
                    rot_matrix[4] = cos(theta) + pow(w[1],2.0) * (1.0 - cos(theta));
                    rot_matrix[5] = w[1] * w[2] * (1.0 - cos(theta)) + w[0] * sin(theta);
                    rot_matrix[6] = w[2] * w[0] * (1.0 - cos(theta)) + w[1] * sin(theta);
                    rot_matrix[7] = w[2] * w[1] * (1.0 - cos(theta)) - w[0] * sin(theta);
                    rot_matrix[8] = cos(theta) + pow(w[2],2.0) * (1.0 - cos(theta));
                }
            }

        public:
            //! Initialize a vector with all 0 elements.
            Vec(void) {
                for (unsigned int i=0; i<dim; ++i) _x[i] = 0;
            };
            //! Initialize a vector of size dim with element values specified by the input array.
            Vec(const double vals[dim]) {
                for (unsigned int i=0; i<dim; ++i) _x[i] = vals[i];
            };
            //! Initialize a vector of dim>=1 as (x,0,...0). If dim < 1 an exception is thrown.
            Vec(const double &x) throw(std::out_of_range) {
                if (dim < 1) throw std::out_of_range("Vec");

                _x[0] = x;
            };
            //! Initialize a vector of dim>=2 as (x,y,0,...0). If dim < 2 an exception is thrown.
            Vec(const double &x, const double &y) throw(std::out_of_range) {
                if (dim < 2) throw std::out_of_range("Vec");

                _x[0] = x;
                _x[1] = y;
            };
            //! Initialize a vector of dim>=3 as (x,y,z,0...0). If dim < 3 an exception is thrown.
            Vec(const double &x, const double &y, const double &z) throw(std::out_of_range) {
                if (dim < 3) throw std::out_of_range("Vec");

                _x[0] = x;
                _x[1] = y;
                _x[2] = z;
            };

            static Vec<dim> nan_vec(void) {
                Vec<dim> new_vec;

                for (unsigned int i=0; i<dim; ++i) new_vec[i] = std::numeric_limits<double>::quiet_NaN();

                return new_vec;
            };

            //! Returns the cross product between this vector and the specified vector.
            //! TODO: Keep this?  If the result is below a specified tolerance, returns null vector
            Vec<dim> cross(const Vec<dim> &vec) const throw(std::domain_error) {
                if (dim != 3) throw std::domain_error("Vec::cross");

                double vals[dim];
                vals[0] = _x[1] * vec._x[2] - _x[2] * vec._x[1];
                vals[1] = _x[2] * vec._x[0] - _x[0] * vec._x[2];
                vals[2] = _x[0] * vec._x[1] - _x[1] * vec._x[0];
                Vec<dim> res(vals);

                if (res.mag() >= CROSS_TOLERANCE) return res;
                else return Vec<dim>();
            };

            //! Calculate the angle (in radians) between this vector and another.
            double vector_angle(const Vec<dim> &vec) const {
                double cos_r = this->dot_product(vec)/(this->mag() * vec.mag());

                if (cos_r > 1.0) cos_r = 1.0;

                if (cos_r < -1.0) cos_r = -1.0;

                return acos(cos_r);
            };

            //! Calculate the distance between this point and another.
            double dist(const Vec<dim> &a) const {
                return (*this-a).mag();
            };

            //! Calculate the magnitude of this vector.
            double mag(void) const {
                double v=0;

                for (unsigned int i=0; i<dim; ++i) v += pow(_x[i], 2);

                return sqrt(v);
            };

            //! Calculate the dot product of this vector with another.
            double dot_product(const Vec<dim> &a) const {
                double v=0;

                for (unsigned int i=0; i<dim; ++i) v += _x[i]*a._x[i];

                return v;
            };

            //! Get the unit vector along this vector.
            Vec<dim> unit_vector(void) const {
                double m = mag(), vals[dim];

                for (unsigned int i=0; i<dim; ++i) vals[i] = _x[i]/m;

                return Vec<dim>(vals);
            };

            //! Add a vector to this vector.
            Vec<dim> &operator+= (const Vec<dim> &a) {
                for (unsigned int i=0; i<dim; ++i) _x[i] += a._x[i];

                return *this;
            };
            //! Add a vector to this vector.
            const Vec<dim> operator+ (const Vec<dim> &a) const {
                Vec<dim> res(*this);
                res += a;
                return res;
            };

            //! Subtract a vector from this vector.
            Vec<dim> &operator-= (const Vec<dim> &a) {
                for (unsigned int i=0; i<dim; ++i) _x[i] -= a._x[i];

                return *this;
            };
            //! Subtract a vector from this vector.
            const Vec<dim> operator- (const Vec<dim> &a) const {
                Vec<dim> res(*this);
                res -= a;
                return res;
            };
            //! Unary negative vector operation.
            const Vec<dim> operator- (void) const {
                Vec<dim> res(*this);
                res *= -1.0;
                return res;
            };

            //! Multiply all elements of this vector by a constant.
            void operator*=(const double &a) {
                for (unsigned int i=0; i<dim; ++i) _x[i] *= a;
            };
            //! Multiply all elements of this vector by a constant.
            const Vec<dim> operator* (const double &a) const {
                Vec<dim> res(*this);
                res *= a;
                return res;
            };

            //! Divide all elements of this vector by a constant.
            void operator/=(const double &a) {
                for (unsigned int i=0; i<dim; ++i) _x[i] /= a;
            };
            //! Divide all elements of this vector by a constant.
            const Vec<dim> operator/ (const double &a) const {
                Vec<dim> res(*this);
                res /= a;
                return res;
            };

            //! Test equality of vectors.
            bool operator==(const Vec<dim> &a) {
                for (unsigned int i=0; i<dim; ++i) if (_x[i] != a._x[i]) return false;

                return true;
            };
            //! Test inequality of vectors.
            bool operator!=(const Vec<dim> &a) {
                return !(*this==a);
            };

            //! Subscript operator to read vector elements.
            double operator[](const unsigned int &d) const throw(std::out_of_range) {
                if (dim <= d) throw std::out_of_range("Vec[]");

                return _x[d];
            };
            //! Subscript operator to read/modify vector elements.
            double &operator[](const unsigned int &d) throw(std::out_of_range) {
                if (dim <= d) throw std::out_of_range("Vec[]");

                return _x[d];
            };

            //! Vector rotation related functions
            Vec<dim> rotate_around_axis(Vec<dim> axis, const double &theta) throw(std::domain_error) {
                if (dim != 3) throw std::domain_error("Vec::rotate_around_axis");

                unsigned int    i, n;
                double          r[dim*dim], outv[dim];
                axis.rotation_matrix(r, theta);

                for (i=0; i<dim; ++i) {
                    outv[i] = 0;

                    for (n=0; n<dim; ++n) outv[i] += r[i*dim+n]*_x[n];
                }

                return Vec<dim>(outv);
            }

            // Return number of bytes used by this object
            unsigned long mem_bytes(void) const {
                return sizeof(double)*dim;
            };

            //! Ordering operator for vectors, only for use with STL.
            bool operator< (const Vec<dim> &a) const {
                for (unsigned int i=0; i<dim; ++i) if (_x[i] < a._x[i]) return true;

                return false;
            };
    };

    std::ostream &operator<<(std::ostream &os, const Vec<2> &pt);
    std::ostream &operator<<(std::ostream &os, const Vec<3> &pt);

    typedef std::vector< double > FloatList;
    typedef std::vector< Vec<3> > VectorList;
    typedef std::set<UIndex> SquareIDSet;

    template <unsigned int ncols>
    class TensorRow {
        private:
            double          _entries[ncols];

        public:
            //! Multiply the elements of this row by a vector.
            double operator* (const Vec<ncols> &vec) const {
                double res;
                unsigned int i;

                for (i=0,res=0; i<ncols; ++i) res += vec[i]*_entries[i];

                return res;
            };
            //! Subscript operator to read vector elements.
            const double operator[](const unsigned int &c) const throw(std::out_of_range) {
                if (ncols <= c) throw std::out_of_range("TensorRow[]");

                return _entries[c];
            };
            //! Subscript operator to read/modify vector elements.
            double &operator[](const unsigned int &c) throw(std::out_of_range) {
                if (ncols <= c) throw std::out_of_range("TensorRow[]");

                return _entries[c];
            };
    };

    template <unsigned int ncols, unsigned int nrows>
    class Tensor {
        private:
            TensorRow<ncols>        _rows[nrows];

        public:
            //! Multiply the elements of this tensor by a vector.
            Vec<nrows> operator* (const Vec<ncols> &vec) const {
                Vec<nrows> res;
                unsigned int i;

                for (i=0; i<nrows; ++i) res[i] = _rows[i]*vec;

                return res;
            };
            //! Subscript operator to read vector elements.
            const TensorRow<ncols> operator[](const unsigned int &r) const throw(std::out_of_range) {
                if (nrows <= r) throw std::out_of_range("Tensor[]");

                return _rows[r];
            };
            //! Subscript operator to read/modify vector elements.
            TensorRow<ncols> &operator[](const unsigned int &r) throw(std::out_of_range) {
                if (nrows <= r) throw std::out_of_range("Tensor[]");

                return _rows[r];
            };
    };

    std::ostream &operator<<(std::ostream &os, const Tensor<3,3> &tensor);
    std::ostream &operator<<(std::ostream &os, const TensorRow<3> &tensor);

    /*!
     A latitude/longitude and depth under or on the Earth used in calculations.
     Latitude must be in [-90, 90] degrees and longitude must be within [-180, 180] degrees.
     Altitude indicates meters above ground (negative is underground).
     */
    class LatLonDepth {
        private:
            double      _lat, _lon, _altitude;

        public:
            //! Default constructor - sets latitude, longitude and altitude to be 0.
            LatLonDepth(void) : _lat(std::numeric_limits<double>::quiet_NaN()), _lon(std::numeric_limits<double>::quiet_NaN()), _altitude(std::numeric_limits<double>::quiet_NaN()) {};
            //! Constructor with specified latitude and longitude. Altitude defaults to 0 unless specified.
            LatLonDepth(const double &lat, const double &lon, const double &altitude=0) throw(std::invalid_argument) : _lat(lat), _lon(lon), _altitude(altitude) {
                if (fabs(lat)>90) throw std::invalid_argument("LatLonDepth::lat must be in [-90,90].");

                if (fabs(lon)>180) throw std::invalid_argument("LatLonDepth::lon must be in [-180,180].");
            }

            //! Get the latitude in degrees of this point.
            double lat(void) const {
                return _lat;
            };
            //! Get the longitude in degrees of this point.
            double lon(void) const {
                return _lon;
            };
            //! Get the altitude in meters of this point.
            double altitude(void) const {
                return _altitude;
            };

            //! Set the latitude in degrees of this point.
            void set_lat(const double &lat) throw(std::invalid_argument) {
                if (fabs(lat)>90) throw std::invalid_argument("LatLonDepth::lat must be in [-90,90].");

                _lat = lat;
            };
            //! Set the longitude in degrees of this point.
            void set_lon(const double &lon) throw(std::invalid_argument) {
                if (fabs(lon)>180) throw std::invalid_argument("LatLonDepth::lon must be in [-180,180].");

                _lon = lon;
            };
            //! Set the altitude in meters of this point.
            void set_altitude(const double &altitude) throw(std::invalid_argument) {
                _altitude = altitude;
            };

            // Return number of bytes used by this object
            unsigned long mem_bytes(void) const {
                return sizeof(double)*3;
            };

            //! Tests point equality (identical latitude, longitude and altitude).
            bool operator==(LatLonDepth &pt) const {
                return (_lat==pt._lat && _lon==pt._lon && _altitude==pt._altitude);
            };
            //! Tests point inequality.
            bool operator!=(LatLonDepth &pt) const {
                return (!(*this == pt));
            };
    };

    std::ostream &operator<<(std::ostream &os, const LatLonDepth &pt);

    /*!
     \brief A class to perform conversions between units.

     This class isn't needed because the conversions are particularly complicated,
     but rather to remind the user what units are being converted to what.
     The class also performs longitude/latitude conversion to Cartesian coordinates.
     */
    class Conversion {
        private:
            //! Base coordinate for conversion.
            LatLonDepth         _base;

            void dist_vincenty(double &distance, double &start_azimuth, double &end_azimuth, const LatLonDepth &p1, const LatLonDepth &p2) const;

        public:
            //! Create a Conversion class with a base LatLonDepth of (0,0,0).
            Conversion(void) : _base() {};
            //! Create a Conversion class with the specified base LatLonDepth point.
            Conversion(const LatLonDepth &base) : _base(base) {};

            //! Set the base LatLonDepth used to convert to Cartesian coordinates.
            void set_base_lat_lon_depth(const LatLonDepth &new_base) {
                _base = new_base;
            };
            //! Get the base LatLonDepth.
            LatLonDepth get_base_lat_lon_depth(void) const {
                return _base;
            };

            //! Convert the specified point to a Cartesian coordinate using the base as (0,0,0).
            Vec<3> convert2xyz(const LatLonDepth &in_pt) const;

            //! Convert the specified Cartesian coordinate to latitude/longitude using the base as (0,0,0).
            LatLonDepth convert2LatLon(const Vec<3> &in_pt) const;

            //! Convert degrees to radians.
            static double deg2rad(const double &degrees) {
                return degrees*M_PI/180;
            };

            //! Convert radians to degrees.
            static double rad2deg(const double &radians) {
                return 180*radians/M_PI;
            };

            //! Convert years to seconds.
            static double year2sec(const double &years) {
                return years*365.25*24*60*60;
            };

            //! Convert seconds to years.
            static double sec2year(const double &seconds) {
                return seconds/(60.0*60.0*24.0*365.25);
            };

            //! Convert meters per second into centimeters per year.
            static double m_per_sec2cm_per_yr(const double &mps) {
                return year2sec(mps)*100;
            };

            //! Convert centimeters per year into meters per second.
            static double cm_per_yr2m_per_sec(const double &cpy) {
                return sec2year(cpy)/100;
            };

            //! Convert meters into kilometers.
            static double m2km(const double &meters) {
                return meters*1e-3;
            };

            //! Convert kilometers into meters.
            static double km2m(const double &km) {
                return km*1e3;
            };

            //! Convert centimeters into meters.
            static double cm2m(const double &cm) {
                return cm*1e-2;
            };

            //! Convert square kilometers into square meters.
            static double sqkm2sqm(const double &sqkm) {
                return sqkm*1e6;
            };

            //! Convert square meters into square kilometers.
            static double sqm2sqkm(const double &sqm) {
                return sqm*1e-6;
            };

            //! Convert pascals to bars.
            static double pascal2bar(const double &pascals) {
                return pascals*1e-5;
            };

            //! Convert bars to pascals.
            static double bar2pascal(const double &bars) {
                return bars*1e5;
            };

            VectorList convertArray2xyz(const FloatList &lats, const FloatList &lons) const;
    };


