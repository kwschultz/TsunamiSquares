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
    tsunamisquares::World                  this_world;
    tsunamisquares::Vec<2>                      accel, velo, loc; //auto-init to (0,0)
    tsunamisquares::SquareIDSet::const_iterator it;
    tsunamisquares::SquareIDSet                 ids;
    std::ofstream                               out_file;
    const std::string       file_name = "accel_bump_900_dt10_diffusion.txt";
    // Diffusion constant (fit to a reasonable looking sim)
    double D = 140616.45;
    
    this_world.clear();
    //this_world.read_file_ascii("test_file.txt");
    std::cout << "Reading...  Pacific_900.txt" << std::endl;
    this_world.read_bathymetry("Pacific_900.txt");
    this_world.info();

    ids = this_world.getSquareIDs();

    // Flatten the bathymetry
    double new_depth = -1000.0;
    std::cout << "Flattening the bottom...";
    this_world.flattenBottom(new_depth);
    
    // Put water into squares to bring water level up to sealevel.
    std::cout << "Filling with water...";
    this_world.fillToSeaLevel();
    
    // Initial conditions
    tsunamisquares::UIndex bot_right = (int)(this_world.num_squares()*0.5 + 0.5*sqrt(this_world.num_squares()));
    tsunamisquares::UIndex bot_left  = bot_right-1;
    tsunamisquares::UIndex top_left  = bot_left+(int)sqrt(this_world.num_squares());
    tsunamisquares::UIndex top_right = bot_right+(int)sqrt(this_world.num_squares());
//    // TODO: Save num_lons and num_lats in the world object
    std::cout << "Deforming the bottom... " << std::endl;
    this_world.deformBottom(bot_left,100.0);
    this_world.deformBottom(top_left,100.0);
    this_world.deformBottom(top_right,100.0);
    this_world.deformBottom(bot_right,100.0);
//    this_world.setSquareHeight(0,1001.0);
//    double Lx = this_world.square(0).Lx();
//    double Ly = this_world.square(0).Ly();
//    this_world.setSquareVelocity(0, tsunamisquares::Vec<2>(Lx/2, -Ly/2));
    
    //this_world.write_file_kml("Pacific_36.kml");

    // -------- Prepare a run to write to file ----------------------               
    double dt = (double) (int) this_world.square(0).Lx()*this_world.square(0).Ly()/(2*D); //seconds    
    int N_steps = 20; //number of time steps
    int current_step = 0;
    int update_step = 1;
    int save_step = 1;
    double max_time = N_steps*dt;
    double time = 0.0;
    ids = this_world.getSquareIDs();
    
    // Open the output file
    out_file.open(file_name.c_str());
    // Write the header
    out_file << "# time \t lon \t lat \t height \n";
    std::cout << "Moving squares..";
    while (time < max_time) {
        // If this is a writing step, print status
        if (current_step%update_step == 0) {
            std::cout << ".." << (100.0*current_step)/N_steps << "%..";
            std::cout << std::flush;
        }
    
        // Write the current state to file
        if (current_step%save_step == 0) {
            for (it=ids.begin(); it!=ids.end(); ++it){
                this_world.write_square_ascii(out_file, time, *it);
            }
        }
        // Move the squares
        this_world.moveSquares(dt);
        this_world.diffuseSquares(dt);
        time += dt;
        current_step += 1;
    }
    out_file.close();
    std::cout << std::endl << "Results written to " << file_name << std::endl;
    return 0;
}

