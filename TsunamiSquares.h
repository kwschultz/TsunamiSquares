// Copyright (c) 2015 Kasey W. Schultz, Steven N. Ward, Eric M. Heien
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

// Not going to use any compiler flags/directives yet, assuming HDF5_FOUND=True

#define UNDEFINED_ELEMENT_ID    UINT_MAX

//#include "hdf5.h"
//#include "hdf5_hl.h"
#include <string.h>
#include <fstream>
#include <map>
#include <stdlib.h>

#include "TsunamiSquaresUtil.h"

namespace tsunamisquares {
    typedef unsigned int UIndex;
    static const UIndex INVALID_INDEX = std::numeric_limits<unsigned int>::max();
    
    class ModelIO {
        private:
            std::string         _comment;

        protected:
            std::string comment(void) const {
                return _comment;
            };
            void set_comment(std::string new_comment) {
                _comment = new_comment;
            };
            std::string next_line(std::istream &in_stream);
            void next_line(std::ostream &out_stream) const;
    };

    struct FieldDesc {
        std::string name;
        std::string details;
    };

    typedef struct FieldDesc FieldDesc;

    // Vertices that make up a Tsunami Square
    struct VertexData {
        UIndex  _id;
        double   _lat, _lon, _alt;
        //unsigned int _is_boundary;
    };

    class Vertex : public ModelIO {
        private:
            VertexData          _data;
            Vec<3> _pos;

        public:
            Vertex(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = std::numeric_limits<double>::quiet_NaN();
                _pos = Vec<3>();
                //_data._is_boundary = 0;
            };
            
            void clear(void);

            VertexData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };

            LatLonDepth lld(void) const {
                return LatLonDepth(_data._lat, _data._lon, _data._alt);
            };
            void set_lld(const LatLonDepth &lld, const LatLonDepth &base) {
                Conversion c(base);
                Vec<3> xyz = c.convert2xyz(lld);
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
                _pos = xyz;
            };

            Vec<3> xyz(void) const {
                return _pos;
            };
            Vec<2> xy(void) const {
                Vec<2> _pos_xy;
                _pos_xy[0] = _pos[0];
                _pos_xy[1] = _pos[1];
                return _pos_xy;
            };
            
            void set_xyz(const Vec<3> &new_xyz, const LatLonDepth &base) {
                Conversion c(base);
                LatLonDepth lld = c.convert2LatLon(new_xyz);
                _pos = new_xyz;
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
            };

            void read_bathymetry(std::istream &in_stream);
    };

    // Squares, the functional members of Tsunami Square
    struct SquareData {
        UIndex              _id;
        //unsigned int        _is_boundary;
        // _vertex for the vertex_id for this square
        UIndex              _vertex;
        Vec<2>              _velocity;
        Vec<2>              _accel;
        double              _height;
        double              _friction;
        double              _density;
        double              _area;
        // updated data are used to keep track of volume/momentum redistribution from other squares
        double               _updated_height;
        Vec<2>               _updated_momentum;
    };

    class Square : public ModelIO {
        private:
            SquareData         _data;

        public:
            Square(void) {
                _data._id = INVALID_INDEX;

                _data._vertex = INVALID_INDEX;
                _data._velocity = _data._accel = _data._updated_momentum = Vec<2>(0.0,0.0);

                //_data._is_boundary = false;
                _data._height = _data._updated_height = _data._area  = std::numeric_limits<double>::quiet_NaN();
                _data._density = 1025.0; // sea water by default
                _data._friction = 0.02;
            };
            
            void clear(void);

            SquareData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };
            
            UIndex vertex(void) const {
                return _data._vertex;
            };
            void set_vertex(const UIndex &ind) {
                _data._vertex = ind;
            };
            
            double height(void) const {
                return _data._height;
            };
            void set_height(const double &new_height) {
                _data._height = new_height;
            };
            
            double updated_height(void) const {
                return _data._updated_height;
            };
            void set_updated_height(const double &new_height) {
                _data._updated_height = new_height;
            };
            
            double density(void) const {
                return _data._density;
            };
            void set_density(const double &new_density) {
                _data._density = new_density;
            };
            
            Vec<2> velocity(void) const {
                return _data._velocity;
            };
            void set_velocity(const Vec<2> &new_velocity) {
                _data._velocity = new_velocity;
            };
            
            Vec<2> updated_momentum(void) const {
                return _data._updated_momentum;
            };
            void set_updated_momentum(const Vec<2> &new_velocity) {
                _data._updated_momentum = new_velocity;
            };
            
            Vec<2> accel(void) const {
                return _data._accel;
            };
            void set_accel(const Vec<2> &new_accel) {
                _data._accel = new_accel;
            };
            
            double friction(void) const {
                return _data._friction;
            };
            void set_friction(const double &new_friction) {
                _data._friction = new_friction;
            };
            
            double area(void) const {
                return _data._area;
            };

            void set_area(const double &new_area) {
                _data._area = new_area;
            };
        
            double length(void) const {
                // Compute the side length of the square
                return sqrt(_data._area);
            }

            double volume(void) const {
                return _data._area*_data._height;
            };
    
    };
            
    typedef std::set<UIndex> SquareIDSet;
    
    
    
    // Class to contain all Squares and Bathymetry 
    class World : public ModelIO {
        private:
            std::map<UIndex, Vertex>   _vertices;
            std::map<UIndex, Square>  _squares;
            LatLonDepth _base;
            double _min_lat, _max_lat, _min_lon, _max_lon;
            
        public:
            Square &new_square(void);
            Vertex &new_vertex(void);
            
            Square &square(const UIndex &ind) throw(std::domain_error);
            Vertex &vertex(const UIndex &ind) throw(std::domain_error);
            
            UIndex next_square_index(void) const {
                if (_squares.size()) return _squares.rbegin()->first+1;
                else return 0;
            };
            UIndex next_vertex_index(void) const {
                if (_vertices.size()) return _vertices.rbegin()->first+1;
                else return 0;
            };
            
            size_t num_squares(void) const;
            size_t num_vertices(void) const;
            
            void insert(Square &new_square);
            void insert(const Vertex &new_vertex);
            
            void clear(void);
            
            LatLonDepth getBase(void) const {
                return _base;
            }
            
            void printSquare(const UIndex square_id);
            void printVertex(const UIndex vertex_id);
            void info(void) const;
            
            int read_bathymetry(const std::string &file_name);
            
            void get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const;
            void reset_base_coord(const LatLonDepth &new_base);
            
            SquareIDSet getSquareIDs(void) const;
            SquareIDSet getVertexIDs(void) const;

            SquareIDSet getNearestIDs(const Vec<2> &location) const;
            SquareIDSet getNeighborIDs(const UIndex &square_id) const;
            
            // ======= Main functions =========
            void fillToSeaLevel(void);
            void moveSquares(const double dt);
            void smoothSquares(void);
            Vec<2> getGradient(const UIndex &square_id) const;
            void updateAcceleration(const UIndex &square_id);
            void deformBottom(const UIndex &square_id, const double &height_change);
            UIndex whichSquare(const Vec<2> &location) const;
            void flattenBottom(const double &depth);
            // ======= Square functions =========
            Vec<2> squareCenter(const UIndex &square_id) const;
            Vec<2> squareLatLon(const UIndex &square_id) const;
            double squareDepth(const UIndex &square_id) const;
            double squareLevel(const UIndex &square_id) const;
            double squareMass(const UIndex &square_id) const;
            double squareVolume(const UIndex &square_id) const;
            Vec<2> squareMomentum(const UIndex &square_id) const;
            void write_square_ascii(std::ostream &out_stream, const double &time, const UIndex &square_id) const;
            // ======= Initial condition setup functions ======
            void setSquareVelocity(const UIndex &square_id, const Vec<2> &new_velo);
            void setSquareAccel(const UIndex &square_id, const Vec<2> &new_accel);
            void setSquareHeight(const UIndex &square_id, const double &new_height);
    };
}