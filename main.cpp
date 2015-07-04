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
    tsunamisquares::ModelWorld                  this_world;
    tsunamisquares::Vec<2>                      accel, velo, loc; //auto-init to (0,0)
    tsunamisquares::SquareIDSet::const_iterator it;
    tsunamisquares::SquareIDSet                 ids;
    std::ofstream                               out_file;
    const std::string                           file_name = "test_out.txt";
    
    this_world.clear();
    this_world.read_file_ascii("test_file.txt");
    this_world.info();
    
    // Put water into squares to bring water level up to sealevel.
    this_world.fillToSeaLevel();

    // Give Square 2 a velocity and larger height
    this_world.setSquareVelocity(3,tsunamisquares::Vec<2>(500,500));
    this_world.setSquareHeight(3,2000.0);
    
    float dt = .1; //seconds
    int N_steps = 30; //number of time steps
    float max_time = N_steps*dt;
    float time = 0.0;
    ids = this_world.getSquareIDs();
    
    // Open the output file
    out_file.open(file_name.c_str());
    // Write the header
    out_file << "# time \t square_x \t square_y \t height \n";
    while (time <= max_time) {
        // Write the current state to file
        for (it=ids.begin(); it!=ids.end(); ++it){
            this_world.square(*it).write_ascii_outfile(out_file, time);
        }
        // Move the squares
        this_world.moveSquares(dt);
        time += dt;
    }
    out_file.close();
    std::cout << "Results written to " << file_name << std::endl;
    return 0;
}

