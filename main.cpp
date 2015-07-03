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
    tsunamisquares::Vec<2>      accel, velo, loc; //auto-init to (0,0)
    std::map<double, tsunamisquares::UIndex> dists;
    std::map<double, tsunamisquares::UIndex>::iterator dit;
    tsunamisquares::SquareIDSet::const_iterator it;
    tsunamisquares::SquareIDSet ids;
    
    // Clear the world
    this_world.clear();

    // Read in a model file
    this_world.read_file_ascii("test_file.txt");
    
    // Grab the new square's data from the World
    this_world.info();
    
    // Put water into squares to bring water level up to sealevel.
    this_world.fillToSeaLevel();
    
    // Look at the squares
    ids = this_world.getSquareIDs();
    for (it=ids.begin(); it!=ids.end(); ++it){
        this_world.printSquare(*it);
    }
    
    
    tsunamisquares::Square sq2 = this_world.square(2);
    sq2.set_velocity(tsunamisquares::Vec<2>(1500,1500));
    
    float dt = 1.0; //seconds
    
    this_world.moveSquare(2, dt);
    
    // Look at the squares
    ids = this_world.getSquareIDs();
    for (it=ids.begin(); it!=ids.end(); ++it){
        this_world.printSquare(*it);
    }
    
    return 0;
}

