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
#include <cassert>

#define assertThrow(COND, ERR_MSG) assert(COND);

// ----------------------------------------------------------------------
// -------------------- Main Functions --------------------------------
// ----------------------------------------------------------------------

// Set the height for all elements equal to the depth of the bathymetry below the center of the square.
// Result is squares with just enough water so that the water sits at sea level.
void tsunamisquares::World::fillToSeaLevel(void) {
    std::map<UIndex, Square>::iterator     it;

    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        // Add water if the  altitude of the Square center is below sea level
        if (squareDepth(it->first) < 0.0) { 
            it->second.set_height(fabs(squareDepth(it->first)));
        } else {
            it->second.set_height(0.0);
        }
        // Also initialize the velocities and accel to zero
        it->second.set_velocity(Vec<2>(0.0,0.0));
        it->second.set_accel(Vec<2>(0.0,0.0));
    }
}


// Diffusion: Remove a volume of water from each square and distribute it to the neighbors.
// Model: area_change = diff_const*time_step
void tsunamisquares::World::diffuseSquares(const double dt) {
    std::map<UIndex, Square>::iterator  it;
    SquareIDSet                         neighborIDs;
    std::map<UIndex, Square>::iterator  nit;
    double                              volume_change, new_level, add_height, height_change;
    Vec<2>                              momentum_change;
    SquareIDSet::iterator               id_it;
    bool debug = false;
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //TODO: Check that I do not diffuse into dry squares (wetting them artificially)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    // Initialize updated_heights and momenta, will use this to store the net height and momentum changes
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        it->second.set_updated_height( it->second.height() );
        it->second.set_updated_momentum( it->second.momentum() );
    }

    // Compute the height changes due to diffusion of water to neighbors
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        if (it->second.height() > 0 && squareLevel(it->first) != 0.0) {
            // Compute the new height after diffusing the water by 1 time step
            new_level = squareLevel(it->first)/(1 + D()*dt/it->second.area());
            volume_change = (it->second.area())*(squareLevel(it->first) - new_level);
            //assertThrow(volume_change >= 0, "Volume change should be positive");
            height_change = new_level - squareLevel(it->first);
            // Transfer the proportional amount of momentum
            momentum_change = (it->second.momentum())*volume_change/(it->second.volume());
            
            if (debug) {
                std::cout << "----> Diffusing Square " << it->second.id() << std::endl;
                std::cout << "volume change: " << volume_change << std::endl;
                std::cout << "old level: " << squareLevel(it->first) << std::endl;
                std::cout << "new level: " << new_level << std::endl;
                std::cout << "-> neighbors " << std::endl;
            }
            
            // For continuity, must self-add 1/4 of the volume change to edges and 1/2 to corners.
            // This also balances the momentum distribution.
            int minLat = squareLatLon(it->first)[0] == min_lat();
            int maxLat = squareLatLon(it->first)[0] == max_lat();
            int minLon = squareLatLon(it->first)[1] == min_lon();
            int maxLon = squareLatLon(it->first)[1] == max_lon();
            int cornerSum = minLat + minLon + maxLon + maxLat;    
            if (cornerSum == 1) {
                height_change += volume_change/( it->second.area()*4.0);
            } else if (cornerSum == 2) {
                height_change += volume_change/( it->second.area()*2.0);
            }
            
            // Add the height change to the updated height
            it->second.set_updated_height(it->second.updated_height() + height_change);
        
            neighborIDs = it->second.get_valid_nearest_neighbors();
            
            for (id_it=neighborIDs.begin(); id_it!=neighborIDs.end(); ++id_it) {
                nit = _squares.find(*id_it);
                // Divide up the diffused volume equally amongst neighbors
                add_height = volume_change/( nit->second.area()*4.0);
                nit->second.set_updated_height( nit->second.updated_height() + add_height);
                nit->second.set_updated_momentum( nit->second.updated_momentum() + momentum_change/4.0);
            }
        }
    }
    
    // Reset the heights and velocities based on the changes
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        it->second.set_height( it->second.updated_height() );
        it->second.set_velocity( it->second.updated_momentum() / it->second.mass());
    }
}


// Move the water from a Square given its current velocity and acceleration.
// Partition the volume and momentum into the neighboring Squares.
void tsunamisquares::World::moveSquares(const double dt) {
    std::map<UIndex, Square>::iterator sit;
    bool debug = false;
    
    // Initialize the updated height and velocity to zero. These are the containers
    // used to keep track of the distributed height/velocity from moving squares.
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        sit->second.set_updated_height(0.0);
        Vec<2> m; m[0] = m[1] = 0.0;
        sit->second.set_updated_momentum(m);
        // Set acceleration based on the current slope of the water surface
        updateAcceleration(sit->first);
        
    }
    
    // Now go through each square and move the water, distribute to neighbors
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        Vec<2> current_velo, current_accel, current_pos, new_pos, new_velo;
        std::map<double, tsunamisquares::UIndex> neighbors;
        std::map<double, tsunamisquares::UIndex>::const_iterator nit;
        
        current_pos = squareCenter(sit->first);
        current_velo = sit->second.velocity();
        current_accel = sit->second.accel();
        
        // Move the square
        new_pos = current_pos + current_velo*dt + current_accel*0.5*dt*dt;
        new_velo = current_velo + current_accel*dt;
        
        // If this square moves, distribute the volume and momentum
        if (new_pos!=current_pos) {
        
            if (debug) {
                std::cout << "---Moving Square " << sit->first << std::endl;
                std::cout << "current pos: " << current_pos << std::endl;
                std::cout << "current velo: " << current_velo << std::endl;
                std::cout << "current accel: " << current_accel << std::endl;
                std::cout << "new pos: " << new_pos << std::endl;
                std::cout << "new velo: " << new_velo << std::endl;
            }
        
            // Find the nearest squares and their distances to the new position
            neighbors = getNearest(new_pos);
            //neighbors = getNearest_from(new_pos, sit->first);
            
            // Init these for renormalizing the fractions
            double fraction_sum = 0.0;
            std::map<UIndex, double> originalFractions, renormFractions;
            std::map<UIndex, double>::iterator frac_it;
            std::vector<double> to_erase; // To store the keys of the neighbors map elements to delete
            std::vector<double>::iterator eit;
            
            // Iterate through neighbors once to compute the fractional area overlap.
            for (nit=neighbors.begin(); nit!=neighbors.end(); ++nit) {
                // This iterator will give us the neighbor square 
                std::map<UIndex, Square>::iterator neighbor_it = _squares.find(nit->second);
                double dx = fabs(new_pos[0] - squareCenter(neighbor_it->first)[0]);
                double dy = fabs(new_pos[1] - squareCenter(neighbor_it->first)[1]);
                double Lx = sit->second.Lx();
                double Ly = sit->second.Ly();
                double this_fraction = (1-dx/Lx)*(1-dy/Ly);
                
                // Add neighbors until we get 4 that have fraction > 0 
                if (this_fraction > 0 && originalFractions.size() < 4) {
                    fraction_sum += this_fraction;
                    originalFractions.insert(std::make_pair(nit->second, this_fraction));
                    if (debug) {
                        std::cout << "--neighbor " << nit->second << std::endl;
                        std::cout << "dx/Lx: " << dx/Lx << std::endl;
                        std::cout << "dy/Ly: " << dy/Ly << std::endl;
                        std::cout << "fraction: " << this_fraction << std::endl;
                    }
                } else {
                    // if the fraction is less than 0, then it is not a valid neighbor
                    to_erase.push_back(nit->first);
                }
            }
            
            // Remove invalid neighbors from the neighbors set
            for (eit=to_erase.begin(); eit!=to_erase.end(); ++eit) {
                if (debug) std::cout << "--> ERASING NEIGHBOR " << neighbors.find(*eit)->second << std::endl << std::flush;
                neighbors.erase(*eit);
            }
            
            if (debug) std::cout << "summed (over " << originalFractions.size() << ") Volume fraction: " << fraction_sum << std::endl;
            
            // Then normalize these fractions to enforce conservation.
            for (frac_it=originalFractions.begin(); frac_it!=originalFractions.end(); ++frac_it) {
                //assertThrow((frac_it->second)/fraction_sum <= 1, "Area fraction must be less than 1.");
                renormFractions.insert(std::make_pair(frac_it->first, (frac_it->second)/fraction_sum));
            }
            
            // Check that the normalized fractions sum exactly to 1
            double renormSum = 0.0;
            for (frac_it=renormFractions.begin(); frac_it!=renormFractions.end(); ++frac_it) {
                renormSum += frac_it->second;
            }
            
            if (debug) {
                std::cout.precision(17);
                std::cout << "1st Renorm: summed fraction: " << std::fixed << renormSum << std::endl;
            }
            
            // Compute height and momentum imparted to neighbors
            for (nit=neighbors.begin(); nit!=neighbors.end(); ++nit) {
                // This iterator will give us the neighbor square 
                std::map<UIndex, Square>::iterator neighbor_it = _squares.find(nit->second);
                // This iterates through the renormalized fractions
                frac_it = renormFractions.find(nit->second);
                double areaFraction = frac_it->second;
                
                // Update the amount of water in the neighboring square (conserve volume)
                double dV = sit->second.volume()*areaFraction;
                double H = neighbor_it->second.updated_height();
                double A_n = neighbor_it->second.area();
                neighbor_it->second.set_updated_height(H + dV/A_n);
                
                // Conserve momentum, update the velocity accordingly (at the end)
                Vec<2> dM = new_velo*areaFraction*(sit->second.mass());
                Vec<2> M  = neighbor_it->second.updated_momentum();
                neighbor_it->second.set_updated_momentum(M+dM);
            }
            
        } else {
            // For those squares that don't move, don't change anything.
            sit->second.set_updated_height(sit->second.height());
            sit->second.set_updated_momentum(sit->second.momentum());
        }
    }
    
    // Loop again over squares to set new velocity and height from accumulated height and momentum
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        sit->second.set_height(sit->second.updated_height());
        sit->second.set_velocity(sit->second.updated_momentum()/(sit->second.mass()));
    }
    
}

void tsunamisquares::World::updateAcceleration(const UIndex &square_id) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    Vec<2> grav_accel, friction_accel, gradient;
    double G = 9.80665; //mean gravitational acceleration at Earth's surface [NIST]
    
    // Only accelerate the water in this square IF there is water in this square
    if (square_it->second.height() != 0.0) {
        // gravitational acceleration due to the slope of the water surface
        gradient = getGradient(square_id);
        grav_accel = gradient*G*(-1.0);
        
        // frictional acceleration from fluid particle interaction
        friction_accel = square_it->second.velocity()*(square_it->second.velocity().mag())*(square_it->second.friction())/(-1.0*(square_it->second.height()));
        
        // Set the acceleration
        square_it->second.set_accel(grav_accel + friction_accel);
    } else {
        square_it->second.set_accel( Vec<2>(0.0, 0.0) );
    }
}

tsunamisquares::Vec<2> tsunamisquares::World::fitPointsToPlane(const SquareIDSet &square_ids) {
    // --------------------------------------------------------------------
    // Based on StackOverflow article:
    // http://stackoverflow.com/questions/1400213/3d-least-squares-plane
    // --------------------------------------------------------------------
    std::vector<double>             x_vals, y_vals, z_vals;
    SquareIDSet::const_iterator                      id_it;
    Vec<2>                                        gradient;
    int                             i, N = square_ids.size();
    Vec<9> A;
    Vec<3> b, x;
    
    for (id_it=square_ids.begin(); id_it!=square_ids.end(); ++id_it) {
        x_vals.push_back(squareCenter(*id_it)[0]);
        y_vals.push_back(squareCenter(*id_it)[1]);
        z_vals.push_back(squareLevel(*id_it));
    }
    
    // Build the b vector and the A matrix.
    // Single index for matrix, array style. A[i][j] = A_vec[i*3 + j],  N_cols=3
    for (i=0; i<N; ++i) {
            
        b[0] += x_vals[i]*z_vals[i];
        b[1] += y_vals[i]*z_vals[i];
        b[2] += z_vals[i];
        
        A[0] += x_vals[i]*x_vals[i];
        A[1] += x_vals[i]*y_vals[i];
        A[2] += x_vals[i];
        A[3] += x_vals[i]*y_vals[i];
        A[4] += y_vals[i]*y_vals[i];
        A[5] += y_vals[i];
        A[6] += x_vals[i];
        A[7] += y_vals[i];
        A[8] += 1.0;
    }
    
//    //std::cout << "\npre x: " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
//    std::cout << "\n\nb: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
//    std::cout << "A: " << A[0] << ", " << A[1] << ", " << A[2] << std::endl;
//    std::cout << "   " << A[3] << ", " << A[4] << ", " << A[5] << std::endl;
//    std::cout << "   " << A[6] << ", " << A[7] << ", " << A[8] << std::endl;
    
    // Matrix solver below is adapted from Virtual Quake
    int     j, k;
    double  v, f, sum;
    int     n = 3;

    for (i=0; i<n; ++i) {
        v = A[i+n*i];

        for (j=i+1; j<n; ++j) {
            f = A[i+n*j]/v;

            for (k=0; k<n; ++k) {
                A[k+n*j] -= f*A[k+n*i];
            }

            b[j] -= f*b[i];
        }
    }

    for (i=n-1; i>=0; --i) {
        sum = b[i];

        for (j=i+1; j<n; ++j) {
            sum -= A[j+n*i]*x[j];
        }

        x[i] = sum/A[i+n*i];
        
    }

    std::cout << "\nsolved for x: " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    
    gradient[0] = x[0];
    gradient[1] = x[1];
    return gradient;

}

//tsunamisquares::Vec<2> tsunamisquares::World::getGradient_planeFit(const UIndex &square_id) const {
void tsunamisquares::World::getGradient_planeFit(const UIndex &square_id) {
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    Vec<2> gradient;
    bool debug = false;
    SquareIDSet square_ids_to_fit;
    
    // Build vector pointers
//    double *A_fit = new double[9];
//    double *b_fit = new double[3];
//    double *x_fit = new double[3];
    
    square_ids_to_fit = square_it->second.get_nearest_neighbors_and_self();
    
    gradient = fitPointsToPlane(square_ids_to_fit);
    
    std::cout << "grabbed gradient = (" << gradient[0] << ", " << gradient[1] << ")" << std::endl;

    
}


tsunamisquares::Vec<2> tsunamisquares::World::getGradient(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    Vec<2> gradient;
    bool debug = false;
    
    // Initialize the 4 points that will be used to approximate the slopes d/dx and d/dy
    // for this square. These are the centers of the neighbor squares.
    Vec<2> center = squareCenter(square_id);   
    
    UIndex leftID   = square_it->second.left();
    UIndex rightID  = square_it->second.right();
    UIndex topID    = square_it->second.top();
    UIndex bottomID = square_it->second.bottom();
    
    // TODO: Better boundary conditions. For now, just set no acceleration along boundary.
    // ALSO: Not accelerating squares whose neighbor is on a boundary.
    if (squareLatLon(square_it->first)[0] == min_lat() || squareLatLon(square_it->first)[0] == max_lat() || squareLatLon(square_it->first)[1] == min_lon() || squareLatLon(square_it->first)[1] == max_lon() || leftID == INVALID_INDEX || rightID == INVALID_INDEX || topID == INVALID_INDEX || bottomID == INVALID_INDEX) {
        gradient = Vec<2>(0.0,0.0);
    } else {
        // Altitude of water level of neighbor squares
        double z_left   = squareLevel(leftID);
        double z_right  = squareLevel(rightID);
        double z_top    = squareLevel(topID);
        double z_bottom = squareLevel(bottomID);
        double z_mid    = squareLevel(square_id);

        // Thickness of water in neighbor squares
        double h_left   = _squares.find(leftID)->second.height();
        double h_right  = _squares.find(rightID)->second.height();
        double h_top    = _squares.find(topID  )->second.height();
        double h_bottom = _squares.find(bottomID)->second.height();
        double h_mid    = square_it->second.height();
        
        // X,Y of neighbor squares
        Vec<2> center_L = squareCenter(leftID);
        Vec<2> center_R = squareCenter(rightID);
        Vec<2> center_T = squareCenter(topID);
        Vec<2> center_B = squareCenter(bottomID);
        
        // ================================================================
        // Gradient = (dz/dx, dz/dy)
        // Handle the cases with dry cells on either left/right/top/bottom.
        // IGNORE cells that are hi and dry
        // ================================================================
        if (h_left == 0.0 && h_right == 0.0 && h_top == 0.0 && h_bottom == 0.0) {
        // Case: No water on any side
        // TODO: Is this check needed?
            gradient[0] = 0.0;
            gradient[1] = 0.0;
        } else  {
            if (h_left > 0.0 && h_right > 0.0 && h_top > 0.0 && h_bottom > 0.0) {
            // Case: No dry neighbors, then do normal gradient
                gradient[0] = (z_right-z_left)/( center_L.dist(center_R) );
                gradient[1] = (z_top-z_bottom)/( center_T.dist(center_B) );
            }
            
            // Case: Hi and dry on the right, water to the left
            if (h_right == 0.0 && z_right >= 0.0 && h_left != 0.0) {
                gradient[0] = (z_mid-z_left)/( center_L.dist(center) );
            } else if (h_left == 0.0 && z_left >= 0.0 && h_right != 0.0) {
            // Case: Hi and dry on the left, water to the right
                gradient[0] = (z_right-z_mid)/( center_R.dist(center) );
            }

            
            // Case: Hi and dry on the top, water on bottom
            if (h_top == 0.0 && z_top >= 0.0 && h_bottom != 0.0) {
                gradient[1] = (z_mid-z_bottom)/( center.dist(center_B) );
            } else if (h_left == 0.0 && z_left >= 0.0 && h_right != 0.0) {
            // Case: Hi and dry on the bottom, water on top
                gradient[1] = (z_top-z_mid)/( center_T.dist(center) );
            }
            
        }
        
        if (debug) {
            std::cout << "square  " << square_id << std::endl;
            std::cout << "d/dx " << gradient[0] << std::endl; 
            std::cout << "d/dy " << gradient[1] << std::endl;
        }
    }
    
    return gradient;
    
}


// Raise/lower the sea floor depth at the square's vertex by an amount "height_change"
void tsunamisquares::World::deformBottom(const UIndex &square_id, const double &height_change) {
    std::map<UIndex, Square>::iterator sit = _squares.find(square_id);
    std::map<UIndex, Vertex>::iterator vit = _vertices.find(sit->second.vertex());
    LatLonDepth new_lld;
    double old_altitude;
    
    new_lld = vit->second.lld();
    old_altitude = new_lld.altitude();
    new_lld.set_altitude(old_altitude + height_change);
    vit->second.set_lld(new_lld, getBase());
}


// Flatten the bottom to be the specified depth
void tsunamisquares::World::flattenBottom(const double &depth) {
    std::map<UIndex, Vertex>::iterator vit;
    LatLonDepth new_lld;
    double newDepth = -fabs(depth);
    
    // Assign the depth for all vertices to be newDepth
    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        new_lld = vit->second.lld();
        new_lld.set_altitude(newDepth);
        vit->second.set_lld(new_lld, getBase());
    }
}
//
//// ----------------------------------------------------------------------
//// -------------------- Utility Functions -------------------------------
//// ----------------------------------------------------------------------

std::map<double, tsunamisquares::UIndex> tsunamisquares::World::getNearest_from(const Vec<2> &location, const UIndex &original_id) const {
    std::map<double, UIndex>                  global_square_dists;
    std::map<double, UIndex>                  local_square_dists;
    SquareIDSet::const_iterator               it, iit;
    std::map<double, UIndex>::const_iterator  dist_it;
    std::map<UIndex, Square>::const_iterator  sit;
    SquareIDSet                               valid_neighbors, valid_neighbors_and_neighbors, neighbors_neighbors; 
    SquareIDSet                               global_neighbors, local_neighbors;
    UIndex                                    neighbor_id;
    
    valid_neighbors = _squares.find(original_id)->second.get_valid_neighbors();
    
    // Iterate over the neighbors of the original position of the square before moving to new location.
    //    Keep track of their IDs, and also grab the IDs of their neighbors. These next-nearest neighbors 
    //    and their neighbors should cover the area around "location" if the time step is small enough.
    Square orig_square = _squares.find(original_id)->second;
    assertThrow(squareCenter(original_id).dist(location)/fmin(orig_square.Lx(),orig_square.Ly()) < 2,"Square moved more than 2 square lengths in one time step!");
    
    //////// Debug -------------------
    std::cout << "\nMoving square " << original_id << " to " << location[0] << "," << location[1] << std::endl;
    
    for (it=valid_neighbors.begin(); it!=valid_neighbors.end(); ++it) {
        neighbor_id = *it;
        SquareIDSet neighbors_neighbors = _squares.find(original_id)->second.get_valid_neighbors();
        // Add the neighbor's neighbors to the cumulative list
        valid_neighbors_and_neighbors.insert(neighbor_id);
        for (iit=neighbors_neighbors.begin(); iit!=neighbors_neighbors.end(); ++iit) {
            valid_neighbors_and_neighbors.insert(*iit);
        }
    }
    
    
    // Compute distance from "location" to the center of each neighboring square.
    for (it=valid_neighbors_and_neighbors.begin(); it!=valid_neighbors_and_neighbors.end(); ++it) {
        double square_dist = squareCenter(*it).dist(location);
        local_square_dists.insert(std::make_pair(square_dist, *it));
    }

    std::cout << "Neighbors nearest:   Square " << local_square_dists.begin()->second << " with distance " << local_square_dists.begin()->first << std::endl;
    std::cout << "Global nearest:   Square " << global_square_dists.begin()->second << " with distance " << global_square_dists.begin()->first << std::endl;

    UIndex local_nearest = local_square_dists.begin()->second;
    UIndex global_nearest = global_square_dists.begin()->second;
    

    std::cout << "Nearest square methods match: " << (local_nearest==global_nearest) << std::endl;
    
    // Compute distance from "location" to the center of each square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        double square_dist = squareCenter(sit->first).dist(location);
        global_square_dists.insert(std::make_pair(square_dist, sit->second.id()));
    }
    
    // Iterate again thru the distance-sorted map, grab the closest squares
//    for (dist_it=global_square_dists.begin(); dist_it!=global_square_dists.end(); ++dist_it) {
//        global_neighbors.insert(std::make_pair(dist_it->first, dist_it->second));
//        if (global_neighbors.size() == valid_neighbors_and_neighbors.size()) break;
//    }
    
    return global_square_dists;
}



// Get the square_id for each of the 4 closest squares to some location = (x,y)
std::map<double, tsunamisquares::UIndex> tsunamisquares::World::getNearest(const Vec<2> &location) const {
    std::map<double, UIndex>                  square_dists;
    std::map<double, UIndex>::const_iterator  it;
    std::map<UIndex, Square>::const_iterator  sit;
    std::map<double, UIndex>                  neighbors;

    // Compute distance from "location" to the center of each square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        double square_dist = squareCenter(sit->first).dist(location);
        square_dists.insert(std::make_pair(square_dist, sit->second.id()));
    }
    
    // Iterate again thru the distance-sorted map, grab the closest squares
    for (it=square_dists.begin(); it!=square_dists.end(); ++it) {
        neighbors.insert(std::make_pair(it->first, it->second));
        if (neighbors.size() == 8) break;
    }
    
    return neighbors;
}

// Get the square_id for each closest square to some location = (x,y)
tsunamisquares::UIndex tsunamisquares::World::whichSquare(const Vec<2> &location) const {
    std::map<double, UIndex>                  square_dists;
    std::map<UIndex, Square>::const_iterator  sit;
    UIndex                               neighbor;

    // Compute distance from "location" to the center of each square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        double square_dist = squareCenter(sit->first).dist(location);
        square_dists.insert(std::make_pair(square_dist, sit->second.id()));
    }
    
    // Return the ID of the nearest square
    return square_dists.begin()->second;
}


void tsunamisquares::World::computeNeighbors(void) {
    std::map<UIndex, Square>::iterator                                  it;
    double                                              this_lat, this_lon; 
    bool                                    minLat, minLon, maxLat, maxLon;
    UIndex                       this_id, left, right, top_right, top_left;
    UIndex                          bottom_left, bottom_right, top, bottom;
    
    // Use the in-place element numbering to find the IDs of the neighboring squares.
    // Must handle the border and corner cases and not include off-model neighbors.
    
    for (it=_squares.begin(); it!=_squares.end(); ++it) {
        this_id     = it->first;
        left        = this_id-1;
        right       = this_id+1;
        top         = this_id-num_lons();
        bottom      = this_id+num_lons();
        top_left    = top-1;
        top_right   = top+1;
        bottom_left = bottom-1;
        bottom_right= bottom+1;
        
        this_lat    = squareLatLon(this_id)[0];
        this_lon    = squareLatLon(this_id)[1];
        minLat      = this_lat == min_lat();
        maxLat      = this_lat == max_lat();
        minLon      = this_lon == min_lon();
        maxLon      = this_lon == max_lon();
        
        // Handle the corner and edge cases
        if (! (maxLat || maxLon || minLon || minLat)) {
            // Interior squares
            it->second.set_right(right);
            it->second.set_left(left);
            it->second.set_top(top);
            it->second.set_bottom(bottom);
            it->second.set_top_left(top_left);
            it->second.set_top_right(top_right);
            it->second.set_bottom_left(bottom_left);
            it->second.set_bottom_right(bottom_right);
        } else if (maxLat && minLon) {
            // Top left (North West) corner
            it->second.set_right(right);
            it->second.set_bottom(bottom);
            it->second.set_bottom_right(bottom_right);
        } else if (maxLat && maxLon) {
            // Top right (North East) corner
            it->second.set_left(left);
            it->second.set_bottom(bottom);
            it->second.set_bottom_left(bottom_left);
        } else if (minLat && maxLon) {
            // Bottom right (South East) corner
            it->second.set_left(left);
            it->second.set_top(top);
            it->second.set_top_left(top_left);
        } else if (minLat && minLon) {
            // Bottom left (South West) corner
            it->second.set_right(right);
            it->second.set_top(top);
            it->second.set_top_right(top_right);
        } else if (minLon) {
            // Left (West) border
            it->second.set_right(right);
            it->second.set_top(top);
            it->second.set_bottom(bottom);
            it->second.set_top_right(top_right);
            it->second.set_bottom_right(bottom_right);
        } else if (maxLat) {
            // Top (North) border
            it->second.set_right(right);
            it->second.set_left(left);
            it->second.set_bottom(bottom);
            it->second.set_bottom_left(bottom_left);
            it->second.set_bottom_right(bottom_right);
        } else if (maxLon) {
            // right (East) border
            it->second.set_left(left);
            it->second.set_top(top);
            it->second.set_bottom(bottom);
            it->second.set_top_left(top_left);
            it->second.set_bottom_left(bottom_left);
        } else if (minLat) {
            // Bottom (South) border
            it->second.set_right(right);
            it->second.set_left(left);
            it->second.set_top(top);
            it->second.set_top_left(top_left);
            it->second.set_top_right(top_right);        
        } else {
            std::cout << "Error, no match to any case! (square " << this_id << ")" << std::endl;
        }

    }
    
    
}

// ----------------------------------------------------------------------
// -------------------- Single Square Functions -------------------------
// ----------------------------------------------------------------------
tsunamisquares::Vec<2> tsunamisquares::World::squareCenter(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    std::map<UIndex, Vertex>::const_iterator  vit = _vertices.find(sit->second.vertex());
    // (x,y) of square center
    return vit->second.xy();
}

double tsunamisquares::World::squareDepth(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    std::map<UIndex, Vertex>::const_iterator  vit = _vertices.find(sit->second.vertex());
    // altitude of the sea floor below this square (negative below sea level)
    return vit->second.xyz()[2];
}

double tsunamisquares::World::squareLevel(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    std::map<UIndex, Vertex>::const_iterator  vit = _vertices.find(sit->second.vertex());
    // altitude of the water surface for this square
    // = altitude of sea floor + height of water
    return (vit->second.xyz()[2])+(sit->second.height());
}

tsunamisquares::Vec<2> tsunamisquares::World::squareLatLon(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator sit = _squares.find(square_id);
    std::map<UIndex, Vertex>::const_iterator  vit = _vertices.find(sit->second.vertex());
    Vec<2> centerLatLon;
    
    centerLatLon[0] = vit->second.lld().lat();
    centerLatLon[1] = vit->second.lld().lon();
    
    return centerLatLon;
}

// ----------------------------------------------------------------------
// -------------------- Functions to set initial conditions  ------------
// ----------------------------------------------------------------------
void tsunamisquares::World::setSquareVelocity(const UIndex &square_id, const Vec<2> &new_velo) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_velocity(new_velo);
}

void tsunamisquares::World::setSquareAccel(const UIndex &square_id, const Vec<2> &new_accel) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_accel(new_accel);
}

void tsunamisquares::World::setSquareHeight(const UIndex &square_id, const double &new_height) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    square_it->second.set_height(new_height);
}

// ----------------------------------------------------------------------
// -------------------- Model Building/Editing --------------------------
// ----------------------------------------------------------------------
tsunamisquares::SquareIDSet tsunamisquares::World::getSquareIDs(void) const {
    SquareIDSet square_id_set;
    std::map<UIndex, Square>::const_iterator  sit;

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        square_id_set.insert(sit->second.id());
    }

    return square_id_set;
}

tsunamisquares::SquareIDSet tsunamisquares::World::getVertexIDs(void) const {
    SquareIDSet vertex_id_set;
    std::map<UIndex, Vertex>::const_iterator  vit;

    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        vertex_id_set.insert(vit->second.id());
    }

    return vertex_id_set;
}

tsunamisquares::Square &tsunamisquares::World::square(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, Square>::iterator it = _squares.find(ind);

    if (it == _squares.end()) throw std::domain_error("tsunamisquares::World::square");
    else return it->second;
}

tsunamisquares::Vertex &tsunamisquares::World::vertex(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, Vertex>::iterator it = _vertices.find(ind);

    if (it == _vertices.end()) throw std::domain_error("tsunamisquares::World::vertex");
    else return it->second;
}

tsunamisquares::Square &tsunamisquares::World::new_square(void) {
    UIndex  max_ind = next_square_index();
    _squares.insert(std::make_pair(max_ind, Square()));
    _squares.find(max_ind)->second.set_id(max_ind);
    return _squares.find(max_ind)->second;
}

tsunamisquares::Vertex &tsunamisquares::World::new_vertex(void) {
    UIndex  max_ind = next_vertex_index();
    _vertices.insert(std::make_pair(max_ind, Vertex()));
    _vertices.find(max_ind)->second.set_id(max_ind);
    return _vertices.find(max_ind)->second;
}

void tsunamisquares::World::clear(void) {
    _squares.clear();
    _vertices.clear();
}

void tsunamisquares::World::reset_base_coord(const LatLonDepth &new_base) {
    std::map<UIndex, Vertex>::iterator         it;

    for (it=_vertices.begin(); it!=_vertices.end(); ++it) {
        it->second.set_lld(it->second.lld(), new_base);
    }

    _base = new_base;
}

void tsunamisquares::World::insert(tsunamisquares::Square &new_square) {
    _squares.insert(std::make_pair(new_square.id(), new_square));
}

void tsunamisquares::World::insert(const tsunamisquares::Vertex &new_vertex) {
    _vertices.insert(std::make_pair(new_vertex.id(), new_vertex));
}

size_t tsunamisquares::World::num_squares(void) const {
    return _squares.size();
}

size_t tsunamisquares::World::num_vertices(void) const {
    return _vertices.size();
}

void tsunamisquares::World::printSquare(const UIndex square_id) {
    Square this_square = square(square_id);

    std::cout << "~~~ Square " << this_square.id() << "~~~" << std::endl;
    std::cout << " vertex " << this_square.vertex() << ": lld " << vertex(this_square.vertex()).lld() << std::endl;
    std::cout << "center: " << squareCenter(square_id) << std::endl;
    std::cout << "density: " << this_square.density() << std::endl;
    std::cout << "area: " << this_square.area() << std::endl;
    if (!isnan(this_square.height())) {
        std::cout << "height: " << this_square.height() << std::endl;
        std::cout << "level: " << squareLevel(square_id) << std::endl;
        std::cout << "volume: " << this_square.volume() << std::endl;
        std::cout << "mass: " << this_square.mass() << std::endl;
        std::cout << "velocity: " << this_square.velocity() << std::endl; 
        std::cout << "accel: " << this_square.accel() << std::endl;    
        std::cout << "momentum: " << this_square.momentum() << std::endl; 
    }
}

void tsunamisquares::World::printVertex(const UIndex vertex_id) {
    Vertex this_vert = vertex(vertex_id);
    std::cout << " ~ Vertex " << this_vert.id() << "~" << std::endl;
    std::cout << "position(xyz): " << this_vert.xyz() << std::endl; 
    std::cout << "position(lld): " << this_vert.lld() << std::endl; 
}

void tsunamisquares::World::info(void) const{
    std::cout << "World: " << this->num_squares() << " squares, " << this->num_vertices() << " vertices." << std::endl;
}

void tsunamisquares::World::get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const {
    std::map<UIndex, Vertex>::const_iterator    it;
    double      min_lat, min_lon, min_alt;
    double      max_lat, max_lon, max_alt;

    min_lat = min_lon = min_alt = DBL_MAX;
    max_lat = max_lon = max_alt = -DBL_MAX;

    for (it=_vertices.begin(); it!=_vertices.end(); ++it) {
        min_lat = fmin(min_lat, it->second.lld().lat());
        max_lat = fmax(max_lat, it->second.lld().lat());
        min_lon = fmin(min_lon, it->second.lld().lon());
        max_lon = fmax(max_lon, it->second.lld().lon());
        min_alt = fmin(min_alt, it->second.lld().altitude());
        max_alt = fmax(max_alt, it->second.lld().altitude());
    }

    if (min_lat == DBL_MAX || min_lon == DBL_MAX || min_alt == DBL_MAX) {
        minimum = LatLonDepth();
    } else {
        minimum = LatLonDepth(min_lat, min_lon, min_alt);
    }

    if (max_lat == -DBL_MAX || max_lon == -DBL_MAX || max_alt == -DBL_MAX) {
        maximum = LatLonDepth();
    } else {
        maximum = LatLonDepth(max_lat, max_lon, max_alt);
    }
}


// ----------------------------------------------------------------------
// ----------------------------- Model File I/O -------------------------
// ----------------------------------------------------------------------
std::string tsunamisquares::ModelIO::next_line(std::istream &in_stream) {
    std::string line = "";
    size_t      pos;

    do {
        std::getline(in_stream, line);
        _comment = "";
        // Cut off any initial whitespace
        pos = line.find_first_not_of(" \t");

        if (pos != std::string::npos) line = line.substr(pos, std::string::npos);

        // Comment consists of hash mark until the end of the line
        pos = line.find("#");

        if (pos != std::string::npos) _comment = line.substr(pos, std::string::npos);

        // Extract the non-comment part of the line
        line = line.substr(0, line.find("#"));

        // If the line is empty, we keep going
        if (line.length() > 0) break;
    } while (in_stream && !in_stream.eof());

    return line;
}

void tsunamisquares::ModelIO::next_line(std::ostream &out_stream) const {
    if (!_comment.empty()) out_stream << " # " << _comment;

    out_stream << "\n";
}


void tsunamisquares::Vertex::read_bathymetry(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));
    
    ss >> _data._lat; 
    ss >> _data._lon;
    ss >> _data._alt;
}

void tsunamisquares::World::write_square_ascii(std::ostream &out_stream, const double &time, const UIndex &square_id) const {
    unsigned int        i;
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    double waterLevel, waterHeight;

    out_stream << time << "\t";

    //
    for (i=0; i<2; ++i) {
        out_stream << squareLatLon(square_id)[i] << "\t\t";
    }

    // Don't write water level for the hi and dry squares until they take on water
    waterLevel  = squareLevel(square_id);
    waterHeight = square_it->second.height();
    if (waterHeight == 0.0 && waterLevel >= 0.0) {
        out_stream << waterHeight << "\t\t";
    } else {
        out_stream << waterLevel << "\t\t";
    }
    
    // Write the altitude of the bottom too
    out_stream << squareDepth(square_id) << "\t\t";

    next_line(out_stream);
}

int tsunamisquares::World::read_bathymetry(const std::string &file_name) {
    std::ifstream   in_file;
    UIndex          i, j, num_squares, num_vertices, num_lats, num_lons;
    LatLonDepth     min_latlon, max_latlon;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    // Read the first line describing the number of sections, etc
    std::stringstream desc_line(next_line(in_file));
    desc_line >> num_lats;
    desc_line >> num_lons;
    _num_latitudes = num_lats;
    _num_longitudes = num_lons;
    
    // Set the number of vertices and squares
    num_vertices = num_lats*num_lons;
    num_squares = num_vertices;

    // Read vertices
    for (i=0; i<num_vertices; ++i) {
        Vertex     new_vert;
        new_vert.read_bathymetry(in_file);
        new_vert.set_id(i);
        _vertices.insert(std::make_pair(new_vert.id(), new_vert));
    }

    in_file.close();
        
    // Reset the internal Cartesian coordinate system
    get_bounds(min_latlon, max_latlon);
    min_latlon.set_altitude(0);
    reset_base_coord(min_latlon);
    // Keep track of Lat/Lon bounds in the World
    _min_lat = min_latlon.lat();
    _min_lon = min_latlon.lon();
    _max_lat = max_latlon.lat();
    _max_lon = max_latlon.lon();
    
    // Manually assign vertex xyz coords for the vertices given the new base lld
    for (i=0; i<num_vertices; ++i) {
        _vertices[i].set_lld(_vertices[i].lld(), getBase());
    }
    
    // Assign the squares
    for (i=0; i<num_squares; ++i) {
        Square     new_square;
        new_square.set_id(i);
        // Assign vertex
        new_square.set_vertex(i);
        // Assign the area from the distance to neighboring vertices.
        // Cases are for the edges of the model.
        double Lx, Ly;
        int row = (int)(i/num_lons);
        int col = (int)(i%num_lons);
        if (col==num_lons-1 && row!=num_lats-1) {
            // right edge, not bottom row
            Lx = (_vertices[i].xy() - _vertices[i-1].xy()).mag();
            Ly = (_vertices[i].xy() - _vertices[i+num_lons].xy()).mag();
        } else if (row==num_lats-1 && col!=num_lons-1) {
            // bottom row not right edge
            Ly = (_vertices[i].xy() - _vertices[i-num_lons].xy()).mag();
            Lx = (_vertices[i].xy() - _vertices[i+1].xy()).mag();
        } else if (row==num_lats-1 && col==num_lons-1) {
            // bottom right corner
            Ly = (_vertices[i].xy() - _vertices[i-num_lons].xy()).mag();
            Lx = (_vertices[i].xy() - _vertices[i-1].xy()).mag();
        } else {
            Lx = (_vertices[i].xy() - _vertices[i+1].xy()).mag();
            Ly = (_vertices[i].xy() - _vertices[i+num_lons].xy()).mag();
        }
        
        new_square.set_Lx(Lx);
        new_square.set_Ly(Ly);
        _squares.insert(std::make_pair(new_square.id(), new_square));
    }
    

    return 0;
}


int tsunamisquares::World::write_file_kml(const std::string &file_name) {
    std::ofstream                             out_file;
    std::map<UIndex, Square>::const_iterator  sit;
    LatLonDepth                               min_bound, max_bound, center;
    Vec<3>                                    min_xyz, max_xyz;
    double                                    dx, dy, range, Lx, Ly;
    unsigned int                              i;
    double                                    depth = 100; //So the squares are off the surface a bit

    out_file.open(file_name.c_str());

    get_bounds(min_bound, max_bound);
    center = LatLonDepth(max_bound.lat() - (max_bound.lat()-min_bound.lat())/2,
                         max_bound.lon() - (max_bound.lon()-min_bound.lon())/2);
    Conversion c(center);
    min_xyz = c.convert2xyz(min_bound);
    max_xyz = c.convert2xyz(max_bound);
    dx = max_xyz[0]-min_xyz[0];
    dy = max_xyz[1]-min_xyz[1];
    range = fmax(dx, dy) * (1.0/tan(c.deg2rad(30)));

    out_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out_file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    out_file << "<Document>\n";
    out_file << "<LookAt>\n";
    out_file << "\t<latitude>" << center.lat() << "</latitude>\n";
    out_file << "\t<longitude>" << center.lon() << "</longitude>\n";
    out_file << "\t<altitude>0</altitude>\n";
    out_file << "\t<range>" << range << "</range>\n";
    out_file << "\t<tilt>0</tilt>\n";
    out_file << "\t<heading>0</heading>\n";
    out_file << "\t<altitudeMode>absolute</altitudeMode>\n";
    out_file << "</LookAt>\n";
    out_file << "<Style id=\"sectionLabel\">\n";
    out_file << "\t<IconStyle>\n";
    out_file << "\t\t<Icon>\n";
    out_file << "\t\t\t<href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>\n";
    out_file << "\t\t</Icon>\n";
    out_file << "\t</IconStyle>\n";
    out_file << "</Style>\n";

    out_file << "<Folder id=\"squares\">\n";

    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        // Compute the lat/lon/depth of the 4 corners of the square
        LatLonDepth         lld[4];
        Vec<2>              v2, centerXY;
        Vec<3>              v3;
        LatLonDepth         centerLLD, base;

        base        = getBase();
        centerLLD   = _vertices[sit->second.vertex()].lld();    
        centerXY    = squareCenter(sit->first);    
        Conversion  c(base);
        Lx          = sit->second.Lx();
        Ly          = sit->second.Ly();
        // Locate the corners in XYZ, then convert to LLD
        v3      = Vec<3>(centerXY[0]-Lx/2.0, centerXY[1]+Ly/2, 0.0); // top left
        lld[0]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]-Lx/2.0, centerXY[1]-Ly/2, 0.0); // bottom left
        lld[1]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+Lx/2.0, centerXY[1]-Ly/2, 0.0); // bottom right
        lld[2]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+Lx/2.0, centerXY[1]+Ly/2, 0.0); // top left
        lld[3]  = c.convert2LatLon(v3);
        
        // Output the KML format polygon for this section
        out_file << "\t\t<Placemark>\n";
        out_file << "\t\t<description>\n";
        out_file << "Square: " << sit->first << "\n";
        out_file << "LLD: " << squareLatLon(sit->first)[0] << "," << squareLatLon(sit->first)[1] << "," << squareDepth(sit->first) << " [m]\n";
        out_file << "XYZ: " << squareCenter(sit->first)[0] << "," << squareCenter(sit->first)[1] << ","   << squareDepth(sit->first) << " [m]\n";
        out_file << "Area: " << sit->second.area()*pow(10,-6) << "[km^2]\n";
        out_file << "Density: " << sit->second.density() << "[kg/m^3]\n";
        out_file << "\t\t</description>\n";
        out_file << "\t\t\t<styleUrl>#baseStyle</styleUrl>\n";
        out_file << "\t\t\t<Polygon>\n";
        out_file << "\t\t\t\t<extrude>0</extrude>\n";
        out_file << "\t\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n";
        out_file << "\t\t\t\t<outerBoundaryIs>\n";
        out_file << "\t\t\t\t\t<LinearRing>\n";
        out_file << "\t\t\t\t\t\t<coordinates>\n";

        for (unsigned int i=0; i<4; ++i) out_file << "\t\t\t\t\t\t\t" << lld[i].lon() << "," << lld[i].lat() << "," << depth << "\n";

        out_file << "\t\t\t\t\t\t</coordinates>\n";
        out_file << "\t\t\t\t\t</LinearRing>\n";
        out_file << "\t\t\t\t</outerBoundaryIs>\n";
        out_file << "\t\t\t</Polygon>\n";
        out_file << "\t\t</Placemark>\n";
    }

    out_file << "</Folder>\n";
    out_file << "</Document>\n";
    out_file << "</kml>\n";

    out_file.close();

    return 0;
}

int tsunamisquares::World::deformFromFile(const std::string &file_name) {
    std::ifstream   in_file;
    UIndex          i, num_points, mappedID;
    double          dz;
    Vec<2>          location;
    std::map<UIndex, Vertex>::iterator vit;
    std::map<UIndex, Square>::iterator sit;
    LatLonDepth     vertex_lld;

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    // Read the first line describing the number of sections, etc
    std::stringstream desc_line(next_line(in_file));
    desc_line >> num_points;

    // Read the points, find nearest square, deform the bottom
    for (i=0; i<num_points; ++i) {
        Vertex     new_vert;
        new_vert.read_bathymetry(in_file);
        new_vert.set_lld(new_vert.lld(), getBase());
        
        // Get location (x,y) for the lat/lon point and get the altitude change
        location = new_vert.xy();
        dz  = new_vert.lld().altitude();
        
        // Find the closest square, grab its vertex
        mappedID = whichSquare(location);
        sit = _squares.find( mappedID );
        vit = _vertices.find( sit->second.vertex() );
        
        // Get the current LLD data for this closest vertex
        vertex_lld = vit->second.lld();
        
        // Update the altitude of the vertex by the amount dz
        vertex_lld.set_altitude( vertex_lld.altitude() + dz);
        
        // Set the new position
        vit->second.set_lld(vertex_lld, getBase());        
        
    }

    in_file.close();

    return 0;
}
