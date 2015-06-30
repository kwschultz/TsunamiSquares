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

#include "hdf5.h"
#include "hdf5_hl.h"
#include <string.h>
#include <map>

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
        size_t      offset;
        hid_t       type;
        size_t      size;
    };

    typedef struct FieldDesc FieldDesc;

    struct VertexData {
        UIndex  _id;
        float   _lat, _lon, _alt;
        unsigned int _is_boundary;
    };

    class Vertex : public ModelIO {
        private:
            VertexData          _data;
            Vec<3> _pos;

        public:
            Vertex(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = std::numeric_limits<float>::quiet_NaN();
                _pos = Vec<3>;
                _data._is_boundary = 0;
            };

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
            void set_xyz(const Vec<3> &new_xyz, const LatLonDepth &base) {
                Conversion c(base);
                LatLonDepth lld = c.convert2LatLon(new_xyz);
                _pos = new_xyz;
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
            };

//            static std::string hdf5_table_name(void) {
//                return "vertices";
//            };
//            
//            static void get_field_descs(std::vector<FieldDesc> &descs);
//
//            void read_data(const VertexData &in_data);
//            void write_data(VertexData &out_data) const;
//
//            void read_ascii(std::istream &in_stream);
//            void write_ascii(std::ostream &out_stream) const;
    };

    struct SquareData {
        UIndex              _id;
        unsigned int        _is_boundary;
        UIndex              _vertices[4];
        float               _velocity[2];
        float               _accel[2];
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
                for (unsigned int i=0; i<2; ++i) _data._velocity[i] = _data._accel[i] = std::numeric_limits<float>::quiet_NaN();

                _data._is_boundary = false;
                _data._height = _data._friction =_data.density = std::numeric_limits<float>::quiet_NaN();
            };

            SquareData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };
            
            bool is_boundary(void) const {
                return _data._is_boundary;
            };
            void set_is_boundary(const bool &is_boundary) {
                _data._is_boundary = is_boundary;
            };
            
            UIndex vertex(const unsigned int &v) const {
                assert(v<4);
                return _data._vertices[v];
            };
            void set_vertex(const unsigned int &v, const UIndex &ind) {
                assert(v<4);
                _data._vertices[v] = ind;
            };
            
            float height(void) const {
                return _data._height;
            };
            void set_height(const float &new_height) {
                _data._height = new_height;
            };
            
            double area(void) const {
                Vec<3> a,b;
                a=_data._vertices[1]-_data._vertices[0];
                b=_data._vertices[2]-_data._vertices[0];
                return a.cross(b).mag();
            };
            
            float volume(void) const {
                return this->area()*this->height();
            };

            float mass(void) const {
                return this->volume()*this->density();
            };

            float density(void) const {
                return _data._density;
            };
            void set_density(const float &new_density) {
                _data._density = new_density;
            };

            //! Calculates the Euclidean distance between the 3D midpoint of this block and another block.
            double center_distance(const Square &other) const {
                return (other.center() - this->center()).mag();
            };
            
            //! Get center point of this block
            Vec<3> center(void) const {
                Vec<3> c;

                for (unsigned int i=0; i<4; ++i) c += _vert[i];

                return c / 4.0;
            };

//            static std::string hdf5_table_name(void) {
//                return "squares";
//            };
//            
//            static void get_field_descs(std::vector<FieldDesc> &descs);
//            void read_data(const SquareData &in_data);
//            void write_data(SquareData &out_data) const;
//
//            void read_ascii(std::istream &in_stream);
//            void write_ascii(std::ostream &out_stream) const;
    };
    
    // Class to contain all Squares and Bathymetry 
    class ModelWorld : public ModelIO {
        private:
            std::map<UIndex, Vertex>   _vertices;
            std::map<UIndex, Square>  _squares;
            LatLonDepth _base;
            double _min_lat, _max_lat, _min_lon, _max_lon;
            
//            void read_square_hdf5(const hid_t &data_file);
//            void read_vertex_hdf5(const hid_t &data_file);
//            void write_square_hdf5(const hid_t &data_file) const;
//            void write_vertex_hdf5(const hid_t &data_file) const;
            
        public:
            Square &new_square(void);
            Vertex &new_vertex(void);
            
            Square &square(const UIndex &ind) throw(std::domain_error);
            Vertex &vertex(const UIndex &ind) throw(std::domain_error);
            
            siterator begin_square(void) {
                return siterator(&_squares, _squares.begin());
            };
            siterator end_square(void) {
                return siterator(&_squares, _squares.end());
            };
            
            UIndex next_square_index(void) const {
                if (_square.size()) return _square.rbegin()->first+1;
                else return 0;
            };
            UIndex next_vertex_index(void) const {
                if (_vertices.size()) return _vertices.rbegin()->first+1;
                else return 0;
            };
            
            size_t num_squares(void) const;
            size_t num_vertices(void) const;
            
            void insert(const Square &new_square);
            void insert(const Vertex &new_vertex);
            
            void clear(void);
            
            void reset_base_coord(const LatLonDepth &new_base);
            
            SquareIDSet getElementIDs(void) const;
            SquareIDSet getVertexIDs(void) const;

            SquareIDSet neighbors(void) const;
            
//            int read_file_ascii(const std::string &file_name);
//            int write_file_ascii(const std::string &file_name) const;            
//            int read_file_hdf5(const std::string &file_name);
//            int write_file_hdf5(const std::string &file_name) const;

    };
    
    // Iterator class for parsing through square objects
    class siterator {
        private:
            std::map<UIndex, Square>              *_map;
            std::map<UIndex, Square>::iterator    _it;
            UIndex                                      _fid;

        public:
            eiterator(void) : _map(NULL), _fid(INVALID_INDEX) {};
            eiterator(std::map<UIndex, Square> *map, std::map<UIndex, Square>::iterator start, const UIndex &fid) : _map(map), _it(start), _fid(fid) {
                if (_map && _fid != INVALID_INDEX) {
                    while (_it != _map->end() && _it->second.section_id() != _fid) {
                        _it++;
                    }
                }
            };
            eiterator &operator=(const eiterator &other) {
                _map = other._map;
                _it = other._it;
                _fid = other._fid;
                return *this;
            };
            bool operator==(const eiterator &other) {
                return (_map == other._map && _it == other._it && _fid == other._fid);
            };
            bool operator!=(const eiterator &other) {
                return (_map != other._map || _it != other._it || _fid != other._fid);
            };
            eiterator &operator++(void) {
                if (_map && _it != _map->end()) {
                    if (_fid == INVALID_INDEX) {
                        _it++;
                    } else {
                        do {
                            _it++;
                        } while (_it != _map->end() && _it->second.section_id() != _fid);
                    }
                }

                return *this;
            };
            Square &operator*(void) {
                return _it->second;
            };
            Square *operator->(void) {
                return (&*(eiterator)*this);
            };
    };
            
    typedef std::set<UIndex> SquareIDSet;
};
    
    
    