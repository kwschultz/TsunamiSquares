// Copyright (c) 2015 Kasey W. Schultz, Steven N. Ward
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

#include "TsunamiLibUtil.h"

namespace tsunamilib {
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

//    struct FieldDesc {
//        std::string name;
//        std::string details;
//        size_t      offset;
//        hid_t       type;
//        size_t      size;
//    };
//
//    typedef struct FieldDesc FieldDesc;

    struct VertexData {
        UIndex  _id;
        float   _lat, _lon, _alt;
        unsigned int _is_boundary;
    };

    class Vertex : public ModelIO {
        private:
            VertexData          _data;
            std::vector<double> _pos;

        public:
            Vertex(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = std::numeric_limits<float>::quiet_NaN();
                _pos = std::vector<double>(3);
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
                std::vector<double>(3) xyz = c.convert2xyz(lld);
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
                _pos = xyz;
            };

            std::vector<double>(3) xyz(void) const {
                return _pos;
            };
            void set_xyz(const std::vector<double>(3) &new_xyz, const LatLonDepth &base) {
                Conversion c(base);
                LatLonDepth lld = c.convert2LatLon(new_xyz);
                _pos = new_xyz;
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
            };

            static std::string hdf5_table_name(void) {
                return "vertices";
            };
            
            static void get_field_descs(std::vector<FieldDesc> &descs);

            void read_data(const VertexData &in_data);
            void write_data(VertexData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    struct SquareData {
        UIndex              _id;
        unsigned int        _is_boundary;
        UIndex              _vertices[4];
        float               _velocity[2];
        float               _accel[2];
        float               _height;
        float               _friction;
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
                _data._height = _data._friction = std::numeric_limits<float>::quiet_NaN();
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


    
    
    
    
    