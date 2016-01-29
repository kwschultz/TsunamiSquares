// Copyright (c) 2015 Kasey W. Schultz
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
#include <time.h>

int main (int argc, char **argv) {
    // Initialize the world (where the squares live), squares and vertices
    tsunamisquares::World                       this_world;
    tsunamisquares::SquareIDSet::const_iterator it;
    tsunamisquares::SquareIDSet                 ids;
    std::ofstream                               out_file;
    clock_t                                     start,end;
    
    
    // -------------------------------------------------------------------------------- //
    ///////////          CONSTANTS (to be moved to sim parameter file)        ////////////    
    // -------------------------------------------------------------------------------- //
    const std::string       out_file_name    = "local/Channel_Islands_bump_flatBottom.txt";
    const std::string       bathy_file       = "local/Channel_Islands.txt";
    //const std::string       kml_file         = "local/Pacific_900.kml";
    const std::string       deformation_file = "local/Channel_Islands_interp_larger_subset_dispField_event1157.txt";
    const std::string       header           = "# time \t lon \t\t lat \t\t water height \t altitude \n";
    
    // Diffusion constant (fit to a reasonable looking sim)
    double D = 140616.45; //140616.45;
    // Flattening the bathymetry to a constant depth (negative for below sea level)
    double new_depth = -100.0;
    // Number of times to move squares
    int N_steps = 100; //number of time steps
    // Updating intervals, etc.
    int current_step = 0;
    int update_step = 1;
    int save_step = 1;
    double time = 0.0;
    int output_num_digits_for_percent = 3;
    
    
    
    // -------------------------------------------------------------------------------- //
    ///////                Simulation Initialization and Loading                   ///////
    // --------------------------------------------------------------------------------//
    start = clock();
    // Read in the bathymetry data
    this_world.clear();
    std::cout << std::endl << "Reading..."   << bathy_file.c_str() << std::endl;
    this_world.read_bathymetry(bathy_file.c_str());
    
    // Gather model information
    this_world.info();
    int num_lats = (int) sqrt(this_world.num_squares());
    int num_lons = num_lats;
    ids = this_world.getSquareIDs();
    double max_time = N_steps*dt;

    // Write KML model
    //std::cout << "Writing KML..."   << kml_file.c_str() << "  ...";
    //this_world.write_file_kml(kml_file.c_str());
    
    // Compute the time step given the diffusion constant D
    double dt = (double) (int) this_world.square(0).Lx()*this_world.square(0).Ly()/(2*D); //seconds

    // Flatten the bottom for simple simulation test cases, do not do this for tsunami simulations
    std::cout << "Flattening the bottom...";
    this_world.flattenBottom(new_depth);
    
    // Put water into squares to bring water level up to sealevel.
    std::cout << "Filling with water..." << std::flush;
    this_world.fillToSeaLevel();
    
    

    // --------------------------------------------------------------------------------//
    //            Sea Floor Deformation and Initial Conditions                         //
    // --------------------------------------------------------------------------------//
    std::cout << "Deforming the bottom... " << std::endl;
    
    // DEFORM VIA FILE
    //this_world.deformFromFile(deformation_file.c_str());

    // Initial conditions
    tsunamisquares::UIndex bot_right = (int)(this_world.num_squares()*0.5 + 0.5*sqrt(this_world.num_squares()));
    tsunamisquares::UIndex bot_left  = bot_right-1;
    tsunamisquares::UIndex top_left  = bot_left+(int)sqrt(this_world.num_squares());
    tsunamisquares::UIndex top_right = bot_right+(int)sqrt(this_world.num_squares());
    ////    // TODO: Save num_lons and num_lats in the world object
    this_world.deformBottom(bot_left,100.0);
    this_world.deformBottom(top_left,100.0);
    this_world.deformBottom(top_right,100.0);
    this_world.deformBottom(bot_right,100.0);



    // --------------------------------------------------------------------------------//
    // --==                         File I/O Preparation                          --== //   
    // --------------------------------------------------------------------------------//            
    out_file.open(out_file_name.c_str());
    // Write the header
    out_file << header.c_str();
    std::cout.precision(output_num_digits_for_percent);
    
    
    
    // --------------------------------------------------------------------------------//
    // --========-           Begin the Simulation; Move the Squares          ----====- //          
    // --------------------------------------------------------------------------------//    
    std::cout << "Moving squares....time_step=" <<dt << "...";
    while (time < max_time) {
        // If this is a writing step, print status
        if (current_step%update_step == 0) {
            std::cout << ".." << (100.0*current_step)/N_steps << "%..";
            std::cout << std::flush;
        }
    
        // Write the current state to file
        if (current_step%save_step == 0) {
            for (it=ids.begin(); it!=ids.end(); ++it) {
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
    
    
    // --------------------------------------------------------------------------------//
    // --========---                    Wrap up and Reporting            ---=======--- //
    // --------------------------------------------------------------------------------//        
    std::cout << std::endl << "Results written to " << out_file_name << std::endl;
    end = clock();
    std::cout.precision(2+output_num_digits_for_percent);
    std::cout << "Total time: " << (float(end)-float(start))/CLOCKS_PER_SEC << " secs." << std::endl << std::endl;
    return 0;
}


//////////// Chalkboard ////////////////

// Grab all squares and print
//    ids = this_world.getSquareIDs();
//    
//    for (it=ids.begin(); it!=ids.end(); ++it) {
//        this_world.printSquare(*it);
//    }


    // Creating up sloping beach over the bottom 5 rows
//    assertThrow(num_lats == num_lons, "lats and lons mismatch");
//    
//    for (unsigned int i=0; i< (int) this_world.num_squares(); ++i) {
//        int row = (int)(i/num_lons);
//        if (row == num_lats-5) this_world.deformBottom(i, 50);
//        if (row == num_lats-4) this_world.deformBottom(i, 75);
//        if (row == num_lats-3) this_world.deformBottom(i, 90);
//        if (row == num_lats-2) this_world.deformBottom(i, 101);
//        if (row == num_lats-1) this_world.deformBottom(i, 110);
//        
//    }

    
//    for (unsigned int i=0; i< (int) this_world.num_squares(); ++i) {
//        // Create a wave coming from the top down, first row
//        int row = (int)(i/num_lons);
//        int col = (int)(i%num_lons);
//        if (row== num_lats-8 && col > 9 && col < num_lons-10) {
//            // Wave is 1m hi
//            this_world.deformBottom(i, 10);
//            //this_world.setSquareVelocity(i, tsunamisquares::Vec<2>(0.0, -this_world.square(0).Ly()/100));
//        }
//    }
    
    // Initial conditions
//    tsunamisquares::UIndex bot_right = (int)(this_world.num_squares()*0.5 + 0.5*sqrt(this_world.num_squares()));
//    tsunamisquares::UIndex bot_left  = bot_right-1;
//    tsunamisquares::UIndex top_left  = bot_left+(int)sqrt(this_world.num_squares());
//    tsunamisquares::UIndex top_right = bot_right+(int)sqrt(this_world.num_squares());
//////    // TODO: Save num_lons and num_lats in the world object
////    std::cout << "Deforming the bottom... " << std::endl;
//    this_world.deformBottom(bot_left,100.0);
//    this_world.deformBottom(top_left,100.0);
//    this_world.deformBottom(top_right,100.0);
//    this_world.deformBottom(bot_right,100.0);
//    
//    // Smooth the initial bump
//    this_world.diffuseSquares(dt);
//    this_world.diffuseSquares(dt);


