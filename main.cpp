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

#include "TsunamiSquares.h"

int main (int argc, char **argv) {
    // Initialize the world (where the squares live), squares and vertices
    tsunamisquares::ModelWorld  this_world;
    tsunamisquares::Square      this_square, that_square;
    tsunamisquares::Vertex      v1, v2, v3, v4;
    tsunamisquares::LatLonDepth base;
    tsunamisquares::Vec<2>      accel, velo;
    
    this_world.clear();
    base = tsunamisquares::LatLonDepth(0.0,0.0,0.0);
    
    //Create vertices for the square, insert them into ModelWorld
    v1.set_id(0);
    v1.set_lld(tsunamisquares::LatLonDepth(0.0,0.0,-10000.0), base);
    this_world.insert(v1);
    
    v2.set_id(1);
    v2.set_lld(tsunamisquares::LatLonDepth(-0.027,0.0,-10000.0), base);
    this_world.insert(v2);
    
    v3.set_id(2);
    v3.set_lld(tsunamisquares::LatLonDepth(-0.027,0.027,-9000.0), base);
    this_world.insert(v3);
    
    v4.set_id(3);
    v4.set_lld(tsunamisquares::LatLonDepth(0.0,0.027,-9000.0), base);
    this_world.insert(v4);
    
    // Initialize a new Square
    this_square.set_id(0);
    this_square.set_is_boundary(true);
    this_square.set_height(10.0);
    this_square.set_density(1000.0);

    velo[0] = 1.0;
    velo[1] = -1.0;
    accel[0] = 1.0;
    accel[1] = 2.0;
    this_square.set_velocity(velo);
    this_square.set_accel(accel);
    
    // Assign vertices to the new Square
    for (unsigned int i=0; i<4; ++i) this_square.set_vertex(i,i);
    
    // Put new Square into the ModelWorld
    this_world.insert(this_square);
    
    // Grab the new square's data from the World
    std::cout << "ModelWorld: " << this_world.num_squares() << " squares, " << this_world.num_vertices() << " vertices." << std::endl;
    std::cout << "~~~~~~~Square " << this_world.square(0).id() << "~~~~~~" << std::endl;
    for (unsigned int i=0; i<4; ++i) {
        std::cout << "  vertex " << this_world.square(0).vertex(i) << ": " << this_world.square(0).vert(i) << std::endl;
    }
    std::cout << "center: " << this_world.square(0).center() << std::endl;
    std::cout << "density: " << this_world.square(0).density() << std::endl;
    std::cout << "height: " << this_world.square(0).height() << std::endl;
    std::cout << "area: " << this_world.square(0).area() << std::endl;
    std::cout << "volume: " << this_world.square(0).volume() << std::endl;
    std::cout << "mass: " << this_world.square(0).mass() << std::endl;
    std::cout << "velocity: (" << this_world.square(0).velocity()[0] << "," << this_world.square(0).velocity()[1] << ")" << std::endl; 
    std::cout << "accel: (" << this_world.square(0).accel()[0] << "," << this_world.square(0).accel()[1] << ")" << std::endl;    
    std::cout << "momentum: (" << this_world.square(0).momentum()[0] << "," << this_world.square(0).momentum()[1] << ")" << std::endl; 
    
    return 0;
}





