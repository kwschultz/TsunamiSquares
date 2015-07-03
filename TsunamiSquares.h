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
        float   _lat, _lon, _alt;
        //unsigned int _is_boundary;
    };

    class Vertex : public ModelIO {
        private:
            VertexData          _data;
            Vec<3> _pos;

        public:
            Vertex(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = std::numeric_limits<float>::quiet_NaN();
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

            
            static void get_field_descs(std::vector<FieldDesc> &descs);

            void read_data(const VertexData &in_data);
            void write_data(VertexData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    // Squares, the functional members of Tsunami Square
    struct SquareData {
        UIndex              _id;
        //unsigned int        _is_boundary;
        // _vertices for the vertex_id's for this square
        UIndex              _vertices[4];
        // _verts for the vertex (x,y,z) positions for this square
        Vec<3>              _verts[4];
        Vec<2>              _velocity;
        Vec<2>              _accel;
        float               _height;
        float               _friction;
        float               _density;
    };

    class Square : public ModelIO {
        private:
            SquareData         _data;

        public:
            Square(void) {
                _data._id = INVALID_INDEX;

                for (unsigned int i=0; i<4; ++i) _data._vertices[i] = INVALID_INDEX;
                for (unsigned int i=0; i<4; ++i) _data._verts[i] = Vec<3>();
                _data._velocity = _data._accel = Vec<2>();

                //_data._is_boundary = false;
                _data._height = _data._friction = std::numeric_limits<float>::quiet_NaN();
                _data._density = 1025.0; // sea water by default
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
            
//            bool is_boundary(void) const {
//                return _data._is_boundary;
//            };
//            void set_is_boundary(const bool &is_boundary) {
//                _data._is_boundary = is_boundary;
//            };
            
            UIndex vertex(const unsigned int &v) const {
                assert(v<4);
                return _data._vertices[v];
            };
            void set_vertex(const unsigned int &v, const UIndex &ind) {
                assert(v<4);
                _data._vertices[v] = ind;
            };

            Vec<3> vert(const unsigned int &v) const {
                assert(v<4);
                return _data._verts[v];
            };
            void set_vert(const unsigned int &v, const Vertex &vertex) {
                assert(v<4);
                _data._verts[v] = vertex.xyz();
            };
            
            float height(void) const {
                return _data._height;
            };
            void set_height(const float &new_height) {
                _data._height = new_height;
            };
            
            float density(void) const {
                return _data._density;
            };
            void set_density(const float &new_density) {
                _data._density = new_density;
            };
            
            Vec<2> velocity(void) const {
                return _data._velocity;
            };
            void set_velocity(const Vec<2> &new_velocity) {
                _data._velocity = new_velocity;
            };
            
            Vec<2> accel(void) const {
                return _data._accel;
            };
            void set_accel(const Vec<2> &new_accel) {
                _data._accel = new_accel;
            };
            
            float friction(void) const {
                return _data._friction;
            };
            void set_friction(const float &new_friction) {
                _data._friction = new_friction;
            };
            
            double area(void) const {
                Vec<3> a,b;
                // Compute area from vertex (x,y) position not using altitude
                a[0] = _data._verts[1][0]-_data._verts[0][0];
                a[1] = _data._verts[1][1]-_data._verts[0][1];
                b[0] = _data._verts[2][0]-_data._verts[0][0];
                b[1] = _data._verts[2][1]-_data._verts[0][1];
                a[2]=b[2]=0.0;
                return a.cross(b).mag();
            };
        
            double length(void) const {
                // Compute the side length of the square
                return std::sqrt(this->area());
            }

            float volume(void) const {
                return this->area()*this->height();
            };

            float mass(void) const {
                return this->volume()*this->density();
            };

            Vec<2> momentum(void) const {
                return _data._velocity*this->mass();
            };

            //! Get center point of this square (at sealevel so z=0)
            Vec<2> center(void) const {
                Vec<3> c;
                Vec<2> cent;
                for (unsigned int i=0; i<4; ++i) c += _data._verts[i];
                cent[0] = c[0];
                cent[1] = c[1];
                return cent / 4.0;
            };

            //! Calculates the Euclidean distance between the midpoint of this block and another block.
            double center_distance(const Square &other) const {
                return (other.center() - this->center()).mag();
            };

            //! Calculates the depth (positive) at the center of the square using the z-coords from the vertices
            double center_depth(void) const {
                double depth = 0.0;
                for (unsigned int i=0; i<4; ++i) depth += _data._verts[i][2];
                return fabs(depth)/4.0;
            };
    
            static void get_field_descs(std::vector<FieldDesc> &descs);
            void read_data(const SquareData &in_data);
            void write_data(SquareData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };
            
    typedef std::set<UIndex> SquareIDSet;
    
    
    
    // Class to contain all Squares and Bathymetry 
    class ModelWorld : public ModelIO {
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
            
            int read_file_ascii(const std::string &file_name);
            int write_file_ascii(const std::string &file_name) const; 
            
            void get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const;
            void reset_base_coord(const LatLonDepth &new_base);
            
            SquareIDSet getSquareIDs(void) const;
            SquareIDSet getVertexIDs(void) const;

            std::map<double, UIndex> getNeighborIDs(const Vec<2> &location) const;
            void fillToSeaLevel(void);
            void moveSquare(const UIndex &square_id);
    };
}