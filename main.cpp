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
    tsunamisquares::SquareIDSet ids; 
    std::set<tsunamisquares::UIndex>::const_iterator it;
    
    this_world.clear();

    // Read in a model file
    this_world.read_file_ascii("test_file.txt");
    
    // Grab the new square's data from the World
    this_world.info();
    
    //loc = this_world.vertex(5).xy();
    //ids = this_world.getNeighborIDs(loc);
    
    // Put water into squares to bring water level up to sealevel.
    this_world.fillToSeaLevel();

    ids = this_world.getSquareIDs();
    for (it=ids.begin(); it!=ids.end(); ++it){
        this_world.printSquare(*it);
    }

    // Try to save the ModelWorld to a file
    //this_world.write_file_ascii("test_file.txt");
    
    return 0;
}