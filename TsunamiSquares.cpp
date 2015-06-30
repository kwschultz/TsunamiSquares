// Copyright (c) 2015 Kasey W. Schultz, Eric M. Heien, Steven N. Ward
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

#include "TsunamiSquares.h"
#include "TsunamiSquaresUtil.h"

tsunamisquares::Square &tsunamisquares::ModelWorld::square(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, Square>::iterator it = _squares.find(ind);

    if (it == _squares.end()) throw std::domain_error("tsunamisquares::ModelWorld::square");
    else return it->second;
}

tsunamisquares::ModelVertex &tsunamisquares::ModelWorld::vertex(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, Vertex>::iterator it = _vertices.find(ind);

    if (it == _vertices.end()) throw std::domain_error("tsunamisquares::ModelWorld::vertex");
    else return it->second;
}

//std::string tsunamisquares::ModelIO::next_line(std::istream &in_stream) {
//    std::string line = "";
//    size_t      pos;
//
//    do {
//        std::getline(in_stream, line);
//        _comment = "";
//        // Cut off any initial whitespace
//        pos = line.find_first_not_of(" \t");
//
//        if (pos != std::string::npos) line = line.substr(pos, std::string::npos);
//
//        // Comment consists of hash mark until the end of the line
//        pos = line.find("#");
//
//        if (pos != std::string::npos) _comment = line.substr(pos, std::string::npos);
//
//        // Extract the non-comment part of the line
//        line = line.substr(0, line.find("#"));
//
//        // If the line is empty, we keep going
//        if (line.length() > 0) break;
//    } while (in_stream && !in_stream.eof());
//
//    return line;
//}

//void tsunamisquares::ModelIO::next_line(std::ostream &out_stream) const {
//    if (!_comment.empty()) out_stream << " # " << _comment;
//
//    out_stream << "\n";
//}

//void tsunamisquares::Square::get_field_descs(std::vector<FieldDesc> &descs) {
//    FieldDesc       field_desc;
//
//    field_desc.name = "id";
//    field_desc.details = "Unique ID of the element.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _id);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "section_id";
//    field_desc.details = "ID of the section associated with the element.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _section_id);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "vertex_0";
//    field_desc.details = "ID of vertex 0.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _vertices[0]);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "vertex_1";
//    field_desc.details = "ID of vertex 1.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _vertices[1]);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "vertex_2";
//    field_desc.details = "ID of vertex 2.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _vertices[2]);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "is_quad";
//    field_desc.details = "Whether the vertices constitute 3 points of a triangle (zero) or 3 points of a parallelogram (non-zero).";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _is_quad);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "slip_rate";
//    field_desc.details = "Long term slip rate of element in meters per second.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _slip_rate);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "aseismic";
//    field_desc.details = "Fraction of slip on element that is aseismic.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _aseismic);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "rake";
//    field_desc.details = "Rake angle of element in radians.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _rake);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "lame_mu";
//    field_desc.details = "Lame's parameter describing the shear modulus of the material for this element (Pascals).";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _lame_mu);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "lame_lambda";
//    field_desc.details = "Lame's lambda parameter of the material for this element, in Pascals.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _lame_lambda);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "max_slip";
//    field_desc.details = "Maximum slip distance for this element, in meters.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _max_slip);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "stress_drop";
//    field_desc.details = "Stress drop for this element, in Pascals.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(ElementData, _stress_drop);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//}
//
//void tsunamisquares::Square::read_data(const ElementData &in_data) {
//    memcpy(&_data, &in_data, sizeof(ElementData));
//}
//
//void tsunamisquares::Square::write_data(ElementData &out_data) const {
//    memcpy(&out_data, &_data, sizeof(ElementData));
//}
//
//void tsunamisquares::Square::read_ascii(std::istream &in_stream) {
//    unsigned int        i;
//    std::stringstream   ss(next_line(in_stream));
//
//    ss >> _data._id;
//    ss >> _data._section_id;
//
//    for (i=0; i<3; ++i) ss >> _data._vertices[i];
//
//    ss >> _data._is_quad;
//    ss >> _data._slip_rate;
//    ss >> _data._aseismic;
//    ss >> _data._rake;
//    ss >> _data._lame_mu;
//    ss >> _data._lame_lambda;
//    ss >> _data._max_slip;
//    ss >> _data._stress_drop;
//}
//
//void tsunamisquares::Square::write_ascii(std::ostream &out_stream) const {
//    unsigned int        i;
//
//    out_stream << _data._id << " ";
//    out_stream << _data._section_id << " ";
//
//    for (i=0; i<3; ++i) out_stream << _data._vertices[i] << " ";
//
//    out_stream << _data._is_quad << " ";
//    out_stream << _data._slip_rate << " ";
//    out_stream << _data._aseismic << " ";
//    out_stream << _data._rake << " ";
//    out_stream << _data._lame_mu << " ";
//    out_stream << _data._lame_lambda << " ";
//    out_stream << _data._max_slip << " ";
//    out_stream << _data._stress_drop;
//
//    next_line(out_stream);
//}

tsunamisquares::Square &tsunamisquares::ModelWorld::new_square(void) {
    UIndex  max_ind = next_square_index();
    _squares.insert(std::make_pair(max_ind, Square()));
    _squares.find(max_ind)->second.set_id(max_ind);
    return _squares.find(max_ind)->second;
}

tsunamisquares::ModelVertex &tsunamisquares::ModelWorld::new_vertex(void) {
    UIndex  max_ind = next_vertex_index();
    _vertices.insert(std::make_pair(max_ind, Vertex()));
    _vertices.find(max_ind)->second.set_id(max_ind);
    return _vertices.find(max_ind)->second;
}

//void tsunamisquares::ModelVertex::get_field_descs(std::vector<FieldDesc> &descs) {
//    FieldDesc       field_desc;
//
//    field_desc.name = "id";
//    field_desc.details = "Unique ID of the vertex.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _id);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "latitude";
//    field_desc.details = "Latitude of the vertex.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _lat);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "longitude";
//    field_desc.details = "Longitude of the vertex.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _lon);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "altitude";
//    field_desc.details = "Altitude of the vertex in meters (negative is below ground).";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _alt);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "das";
//    field_desc.details = "Vertex distance along fault strike in meters.";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _das);
//    field_desc.type = H5T_NATIVE_FLOAT;
//    field_desc.size = sizeof(float);
//#endif
//    descs.push_back(field_desc);
//
//    field_desc.name = "is_trace";
//    field_desc.details = "Whether an element in on the fault trace (non-zero) or not (zero).";
//#ifdef HDF5_FOUND
//    field_desc.offset = HOFFSET(VertexData, _is_trace);
//    field_desc.type = H5T_NATIVE_UINT;
//    field_desc.size = sizeof(UIndex);
//#endif
//    descs.push_back(field_desc);
//
//}
//
//void tsunamisquares::ModelVertex::read_data(const VertexData &in_data) {
//    memcpy(&_data, &in_data, sizeof(VertexData));
//}
//
//void tsunamisquares::ModelVertex::write_data(VertexData &out_data) const {
//    memcpy(&out_data, &_data, sizeof(VertexData));
//}
//
//void tsunamisquares::ModelVertex::read_ascii(std::istream &in_stream) {
//    std::stringstream   ss(next_line(in_stream));
//
//    ss >> _data._id;
//    ss >> _data._lat;
//    ss >> _data._lon;
//    ss >> _data._alt;
//    ss >> _data._das;
//    ss >> _data._is_trace;
//}
//
//void tsunamisquares::ModelVertex::write_ascii(std::ostream &out_stream) const {
//    out_stream << _data._id << " ";
//    out_stream << _data._lat << " ";
//    out_stream << _data._lon << " ";
//    out_stream << _data._alt << " ";
//    out_stream << _data._das << " ";
//    out_stream << _data._is_trace;
//    next_line(out_stream);
//}

void tsunamisquares::ModelWorld::clear(void) {
    _sections.clear();
    _elements.clear();
    _vertices.clear();
}

//int tsunamisquares::ModelWorld::read_file_ascii(const std::string &file_name) {
//    std::ifstream       in_file;
//    unsigned int        i, num_sections, num_elements, num_vertices;
//    LatLonDepth         min_latlon, max_latlon;
//
//    // Clear the world first to avoid incorrectly mixing indices
//    clear();
//
//    in_file.open(file_name.c_str());
//
//    if (!in_file.is_open()) return -1;
//
//    // Read the first line describing the number of sections, etc
//    std::stringstream desc_line(next_line(in_file));
//    desc_line >> num_sections;
//    desc_line >> num_elements;
//    desc_line >> num_vertices;
//
//    // Read sections
//    for (i=0; i<num_sections; ++i) {
//        ModelSection     new_section;
//        new_section.read_ascii(in_file);
//        _sections.insert(std::make_pair(new_section.id(), new_section));
//    }
//
//    // Read elements
//    for (i=0; i<num_elements; ++i) {
//        Square     new_elem;
//        new_elem.read_ascii(in_file);
//        _elements.insert(std::make_pair(new_elem.id(), new_elem));
//    }
//
//    // Read vertices
//    for (i=0; i<num_vertices; ++i) {
//        ModelVertex     new_vert;
//        new_vert.read_ascii(in_file);
//        _vertices.insert(std::make_pair(new_vert.id(), new_vert));
//    }
//
//    in_file.close();
//
//    // Reset the internal Cartesian coordinate system
//    get_bounds(min_latlon, max_latlon);
//    min_latlon.set_altitude(0);
//    reset_base_coord(min_latlon);
//    // Keep track of Lat/Lon bounds in the ModelWorld
//    _min_lat = min_latlon.lat();
//    _min_lon = min_latlon.lon();
//    _max_lat = max_latlon.lat();
//    _max_lon = max_latlon.lon();
//
//    return 0;
//}


void tsunamisquares::ModelWorld::reset_base_coord(const LatLonDepth &new_base) {
    std::map<UIndex, ModelVertex>::iterator         it;

    for (it=_vertices.begin(); it!=_vertices.end(); ++it) {
        it->second.set_lld(it->second.lld(), new_base);
    }

    _base = new_base;
}


//int tsunamisquares::ModelWorld::write_file_ascii(const std::string &file_name) const {
//    std::ofstream                                   out_file;
//    std::vector<FieldDesc>                          descs;
//    std::vector<FieldDesc>::iterator                dit;
//    std::map<UIndex, ModelVertex>::const_iterator   vit;
//    std::map<UIndex, Square>::const_iterator  eit;
//    std::map<UIndex, ModelSection>::const_iterator  fit;
//
//    out_file.open(file_name.c_str());
//    out_file << "# Number of sections\n";
//    out_file << "# Number of elements\n";
//    out_file << "# Number of vertices\n";
//    out_file << _sections.size() << " " << _elements.size() << " " << _vertices.size();
//    next_line(out_file);
//
//    // Write section header
//    descs.clear();
//    ModelSection::get_field_descs(descs);
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << "# " << dit->name << ": " << dit->details << "\n";
//    }
//
//    out_file << "# ";
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << dit->name << " ";
//    }
//
//    out_file << "\n";
//
//    // Write sections
//    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
//        fit->second.write_ascii(out_file);
//    }
//
//    // Write element header
//    descs.clear();
//    Square::get_field_descs(descs);
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << "# " << dit->name << ": " << dit->details << "\n";
//    }
//
//    out_file << "# ";
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << dit->name << " ";
//    }
//
//    out_file << "\n";
//
//    // Write elements
//    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
//        eit->second.write_ascii(out_file);
//    }
//
//    // Write vertex header
//    descs.clear();
//    ModelVertex::get_field_descs(descs);
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << "# " << dit->name << ": " << dit->details << "\n";
//    }
//
//    out_file << "# ";
//
//    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
//        out_file << dit->name << " ";
//    }
//
//    out_file << "\n";
//
//    // Write vertices
//    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
//        vit->second.write_ascii(out_file);
//    }
//
//    out_file.close();
//
//    return 0;
//}
//
//int tsunamisquares::ModelWorld::read_file_hdf5(const std::string &file_name) {
//#ifdef HDF5_FOUND
//    hid_t       plist_id, data_file;
//    herr_t      res;
//    LatLonDepth min_latlon, max_latlon;
//
//    // Clear the world first to avoid incorrectly mixing indices
//    clear();
//
//    if (!H5Fis_hdf5(file_name.c_str())) return -1;
//
//    plist_id = H5Pcreate(H5P_FILE_ACCESS);
//
//    if (plist_id < 0) exit(-1);
//
//    data_file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, plist_id);
//
//    if (data_file < 0) exit(-1);
//
//    read_section_hdf5(data_file);
//    read_element_hdf5(data_file);
//    read_vertex_hdf5(data_file);
//
//    // Release HDF5 handles
//    res = H5Pclose(plist_id);
//
//    if (res < 0) exit(-1);
//
//    res = H5Fclose(data_file);
//
//    if (res < 0) exit(-1);
//
//    // Reset the internal Cartesian coordinate system
//    get_bounds(min_latlon, max_latlon);
//    min_latlon.set_altitude(0);
//    reset_base_coord(min_latlon);
//    // Keep track of Lat/Lon bounds in the ModelWorld
//    _min_lat = min_latlon.lat();
//    _min_lon = min_latlon.lon();
//    _max_lat = max_latlon.lat();
//    _max_lon = max_latlon.lon();
//#else
//    // TODO: Error out
//#endif
//    return 0;
//}
//
//#ifdef HDF5_FOUND
//void tsunamisquares::ModelWorld::read_section_hdf5(const hid_t &data_file) {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, ModelSection>::const_iterator  fit;
//    hsize_t                     num_fields, num_sections;
//    unsigned int                i;
//    SectionData                 *section_data;
//    size_t                      *field_offsets;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    descs.clear();
//    ModelSection::get_field_descs(descs);
//    num_fields = descs.size();
//    field_offsets = new size_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_offsets[i] = descs[i].offset;
//        field_sizes[i] = descs[i].size;
//    }
//
//    res = H5TBget_table_info(data_file, ModelSection::hdf5_table_name().c_str(), &num_fields, &num_sections);
//
//    if (res < 0) exit(-1);
//
//    // TODO: check that num_fields matches the descs
//    //
//    section_data = new SectionData[num_sections];
//    res = H5TBread_records(data_file, ModelSection::hdf5_table_name().c_str(), 0, num_sections, sizeof(SectionData), field_offsets, field_sizes, section_data);
//
//    if (res < 0) exit(-1);
//
//    // Read section data into the World
//    for (i=0; i<num_sections; ++i) {
//        ModelSection  new_section;
//        new_section.read_data(section_data[i]);
//        _sections.insert(std::make_pair(new_section.id(), new_section));
//    }
//
//    // Free memory for HDF5 related data
//    // yoder: ... and use delete [] for arrays...
//    delete [] section_data;
//    delete [] field_offsets;
//    delete [] field_sizes;
//
//}
//
//void tsunamisquares::ModelWorld::read_element_hdf5(const hid_t &data_file) {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, Square>::const_iterator  fit;
//    hsize_t                     num_fields, num_elements;
//    unsigned int                i;
//    ElementData                 *element_data;
//    size_t                      *field_offsets;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    descs.clear();
//    Square::get_field_descs(descs);
//    num_fields = descs.size();
//    field_offsets = new size_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_offsets[i] = descs[i].offset;
//        field_sizes[i] = descs[i].size;
//    }
//
//    res = H5TBget_table_info(data_file, Square::hdf5_table_name().c_str(), &num_fields, &num_elements);
//
//    if (res < 0) exit(-1);
//
//    // TODO: check that num_fields matches the descs
//
//    element_data = new ElementData[num_elements];
//    res = H5TBread_records(data_file, Square::hdf5_table_name().c_str(), 0, num_elements, sizeof(ElementData), field_offsets, field_sizes, element_data);
//
//    if (res < 0) exit(-1);
//
//    // Read element data into the World
//    for (i=0; i<num_elements; ++i) {
//        Square  new_element;
//        new_element.read_data(element_data[i]);
//        _elements.insert(std::make_pair(new_element.id(), new_element));
//    }
//
//    // Free memory for HDF5 related data
//    delete [] element_data;
//    delete [] field_offsets;
//    delete [] field_sizes;
//}
//
//void tsunamisquares::ModelWorld::read_vertex_hdf5(const hid_t &data_file) {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, ModelVertex>::const_iterator  fit;
//    hsize_t                     num_fields, num_vertices;
//    unsigned int                i;
//    VertexData                  *vertex_data;
//    size_t                      *field_offsets;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    descs.clear();
//    ModelVertex::get_field_descs(descs);
//    num_fields = descs.size();
//    field_offsets = new size_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_offsets[i] = descs[i].offset;
//        field_sizes[i] = descs[i].size;
//    }
//
//    res = H5TBget_table_info(data_file, ModelVertex::hdf5_table_name().c_str(), &num_fields, &num_vertices);
//
//    if (res < 0) exit(-1);
//
//    // TODO: check that num_fields matches the descs
//
//    vertex_data = new VertexData[num_vertices];
//    res = H5TBread_records(data_file, ModelVertex::hdf5_table_name().c_str(), 0, num_vertices, sizeof(VertexData), field_offsets, field_sizes, vertex_data);
//
//    if (res < 0) exit(-1);
//
//    // Read vertex data into the World
//    for (i=0; i<num_vertices; ++i) {
//        ModelVertex  new_vertex;
//        new_vertex.read_data(vertex_data[i]);
//        _vertices.insert(std::make_pair(new_vertex.id(), new_vertex));
//    }
//
//    // Free memory for HDF5 related data
//    // yoder: ... and use delete [] for vector/array types...
//    delete [] vertex_data;
//    delete [] field_offsets;
//    delete [] field_sizes;
//}
//
//void tsunamisquares::ModelWorld::write_section_hdf5(const hid_t &data_file) const {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, ModelSection>::const_iterator  fit;
//    size_t                      num_fields, num_sections;
//    unsigned int                i;
//    SectionData                 blank_section, *section_data;
//    char                        **field_names, **field_details;
//    size_t                      *field_offsets;
//    hid_t                       *field_types;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    // Set up the section table definition
//    descs.clear();
//    ModelSection::get_field_descs(descs);
//    num_fields = descs.size();
//    num_sections = _sections.size();
//    field_names = new char *[num_fields];
//    field_details = new char *[num_fields];
//    field_offsets = new size_t[num_fields];
//    field_types = new hid_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_names[i] = new char[descs[i].name.length()+1];
//        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
//        field_names[i][descs[i].name.length()] = '\0';
//        field_details[i] = new char[descs[i].details.length()+1];
//        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
//        field_details[i][descs[i].details.length()] = '\0';
//        field_offsets[i] = descs[i].offset;
//        field_types[i] = descs[i].type;
//        field_sizes[i] = descs[i].size;
//    }
//
//    // TODO: factor this out?
//    blank_section = ModelSection().data();
//
//    // Fill in the data for the sections
//    section_data = new SectionData[num_sections];
//
//    for (i=0,fit=_sections.begin(); fit!=_sections.end(); ++i,++fit) {
//        fit->second.write_data(section_data[i]);
//    }
//
//    // Create the section table
//    res = H5TBmake_table("Fault Sections",
//                         data_file,
//                         ModelSection::hdf5_table_name().c_str(),
//                         num_fields,
//                         num_sections,
//                         sizeof(SectionData),
//                         (const char **)field_names,
//                         field_offsets,
//                         field_types,
//                         num_sections,
//                         &blank_section,
//                         0,
//                         section_data);
//
//    if (res < 0) exit(-1);
//
//    // Add the details of each field as an attribute
//    for (i=0; i<num_fields; ++i) {
//        std::stringstream   ss;
//        ss << "FIELD_" << i << "_DETAILS";
//        res = H5LTset_attribute_string(data_file,
//                                       ModelSection::hdf5_table_name().c_str(),
//                                       ss.str().c_str(),
//                                       field_details[i]);
//
//        if (res < 0) exit(-1);
//    }
//
//    // Free memory for HDF5 related data
//    delete section_data;
//
//    for (i=0; i<num_fields; ++i) delete field_names[i];
//
//    delete field_names;
//
//    for (i=0; i<num_fields; ++i) delete field_details[i];
//
//    delete field_details;
//    delete field_offsets;
//    delete field_types;
//    delete field_sizes;
//}
//
//void tsunamisquares::ModelWorld::write_element_hdf5(const hid_t &data_file) const {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, Square>::const_iterator  eit;
//    size_t                      num_fields, num_elements;
//    unsigned int                i;
//    ElementData                 blank_element, *element_data;
//    char                        **field_names, **field_details;
//    size_t                      *field_offsets;
//    hid_t                       *field_types;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    // Set up the element table definition
//    descs.clear();
//    Square::get_field_descs(descs);
//    num_fields = descs.size();
//    num_elements = _elements.size();
//    field_names = new char *[num_fields];
//    field_details = new char *[num_fields];
//    field_offsets = new size_t[num_fields];
//    field_types = new hid_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_names[i] = new char[descs[i].name.length()+1];
//        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
//        field_names[i][descs[i].name.length()] = '\0';
//        field_details[i] = new char[descs[i].details.length()+1];
//        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
//        field_details[i][descs[i].details.length()] = '\0';
//        field_offsets[i] = descs[i].offset;
//        field_types[i] = descs[i].type;
//        field_sizes[i] = descs[i].size;
//    }
//
//    // Get a blank element for table filling
//    blank_element = Square().data();
//
//    // Fill in the data for the elements
//    element_data = new ElementData[num_elements];
//
//    for (i=0,eit=_elements.begin(); eit!=_elements.end(); ++i,++eit) {
//        eit->second.write_data(element_data[i]);
//    }
//
//    // Create the elements table
//    res = H5TBmake_table("Elements",
//                         data_file,
//                         Square::hdf5_table_name().c_str(),
//                         num_fields,
//                         num_elements,
//                         sizeof(ElementData),
//                         (const char **)field_names,
//                         field_offsets,
//                         field_types,
//                         num_elements,
//                         &blank_element,
//                         0,
//                         element_data);
//
//    if (res < 0) exit(-1);
//
//    // Add the details of each field as an attribute
//    for (i=0; i<num_fields; ++i) {
//        std::stringstream   ss;
//        ss << "FIELD_" << i << "_DETAILS";
//        res = H5LTset_attribute_string(data_file,
//                                       Square::hdf5_table_name().c_str(),
//                                       ss.str().c_str(),
//                                       field_details[i]);
//
//        if (res < 0) exit(-1);
//    }
//
//    // Free memory for HDF5 related data
//    delete element_data;
//
//    for (i=0; i<num_fields; ++i) delete field_names[i];
//
//    delete field_names;
//
//    for (i=0; i<num_fields; ++i) delete field_details[i];
//
//    delete field_details;
//    delete field_offsets;
//    delete field_types;
//    delete field_sizes;
//}
//
//void tsunamisquares::ModelWorld::write_vertex_hdf5(const hid_t &data_file) const {
//    std::vector<FieldDesc>                          descs;
//    std::map<UIndex, ModelVertex>::const_iterator   vit;
//    size_t                      num_fields, num_vertices;
//    unsigned int                i;
//    VertexData                  blank_vertex, *vertex_data;
//    char                        **field_names, **field_details;
//    size_t                      *field_offsets;
//    hid_t                       *field_types;
//    size_t                      *field_sizes;
//    herr_t                      res;
//
//    // Set up the vertex table definition
//    descs.clear();
//    ModelVertex::get_field_descs(descs);
//    num_fields = descs.size();
//    num_vertices = _vertices.size();
//    field_names = new char *[num_fields];
//    field_details = new char *[num_fields];
//    field_offsets = new size_t[num_fields];
//    field_types = new hid_t[num_fields];
//    field_sizes = new size_t[num_fields];
//
//    for (i=0; i<num_fields; ++i) {
//        field_names[i] = new char[descs[i].name.length()+1];
//        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
//        field_names[i][descs[i].name.length()] = '\0';
//        field_details[i] = new char[descs[i].details.length()+1];
//        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
//        field_details[i][descs[i].details.length()] = '\0';
//        field_offsets[i] = descs[i].offset;
//        field_types[i] = descs[i].type;
//        field_sizes[i] = descs[i].size;
//    }
//
//    blank_vertex = ModelVertex().data();
//
//    // Fill in the data for the vertices
//    vertex_data = new VertexData[num_vertices];
//
//    for (i=0,vit=_vertices.begin(); vit!=_vertices.end(); ++i,++vit) {
//        vit->second.write_data(vertex_data[i]);
//    }
//
//    // Create the vertices table
//    res = H5TBmake_table("Vertices",
//                         data_file,
//                         ModelVertex::hdf5_table_name().c_str(),
//                         num_fields,
//                         num_vertices,
//                         sizeof(VertexData),
//                         (const char **)field_names,
//                         field_offsets,
//                         field_types,
//                         num_vertices,
//                         &blank_vertex,
//                         0,
//                         vertex_data);
//
//    if (res < 0) exit(-1);
//
//    // Add the details of each field as an attribute
//    for (i=0; i<num_fields; ++i) {
//        std::stringstream   ss;
//        ss << "FIELD_" << i << "_DETAILS";
//        res = H5LTset_attribute_string(data_file,
//                                       ModelVertex::hdf5_table_name().c_str(),
//                                       ss.str().c_str(),
//                                       field_details[i]);
//
//        if (res < 0) exit(-1);
//    }
//
//    // Free memory for HDF5 related data
//    delete [] vertex_data;
//
//    for (i=0; i<num_fields; ++i) delete [] field_names[i];
//
//    delete [] field_names;
//
//    for (i=0; i<num_fields; ++i) delete [] field_details[i];
//
//    delete [] field_details;
//    delete [] field_offsets;
//    delete [] field_types;
//    delete [] field_sizes;
//}
//#endif
//
//int tsunamisquares::ModelWorld::write_file_hdf5(const std::string &file_name) const {
//#ifdef HDF5_FOUND
//    hid_t       plist_id, data_file;
//    herr_t      res;
//
//    // Create access properties
//    plist_id = H5Pcreate(H5P_FILE_ACCESS);
//
//    if (plist_id < 0) exit(-1);
//
//    // Create the data file, overwriting any old files
//    data_file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//
//    if (data_file < 0) exit(-1);
//
//    write_section_hdf5(data_file);
//    write_element_hdf5(data_file);
//    write_vertex_hdf5(data_file);
//
//    // Release HDF5 handles
//    res = H5Pclose(plist_id);
//
//    if (res < 0) exit(-1);
//
//    res = H5Fclose(data_file);
//
//    if (res < 0) exit(-1);
//
//    return 0;
//#else
//    return 1;
//#endif
//}


void tsunamisquares::ModelWorld::insert(const tsunamisquares::Square &new_square) {
    _squares.insert(std::make_pair(new_square.id(), new_square));
}

void tsunamisquares::ModelWorld::insert(const tsunamisquares::Vertex &new_vertex) {
    _vertices.insert(std::make_pair(new_vertex.id(), new_vertex));
}

size_t tsunamisquares::ModelWorld::num_squares(void) const {
    return return _squares.size();
}

size_t tsunamisquares::ModelWorld::num_vertices(void) const {
    return _vertices.size();
}

tsunamisquares::ElementIDSet tsunamisquares::ModelWorld::getSquareIDs(void) const {
    SquareIDSet square_id_set;
    std::map<UIndex, Square>::const_iterator  sit;

    for (sit=_square.begin(); sit!=_square.end(); ++sit) {
        square_id_set.insert(sit->second.id());
    }

    return square_id_set;
}

tsunamisquares::ElementIDSet tsunamisquares::ModelWorld::getVertexIDs(void) const {
    ElementIDSet vertex_id_set;
    std::map<UIndex, Vertex>::const_iterator  vit;

    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        vertex_id_set.insert(vit->second.id());
    }

    return vertex_id_set;
}