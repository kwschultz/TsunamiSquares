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
        // Add water if the altitude of the Square center is below sea level
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


// Smoothing: Set the volume in each square to the average volume of the neighbors
//void tsunamisquares::World::smoothSquares(void) {
//    std::map<UIndex, Square>::iterator  it;
//    SquareIDSet                         neighborIDs;
//    std::map<UIndex, Square>::iterator  nit;
//    double                              neighbor_volume_avg = 0.0;
//    SquareIDSet::iterator               id_it;
//
//    
//    for (it=_squares.begin(); it!=_squares.end(); ++it) {
//        neighborIDs = getNeighborIDs(it->first);
//        std::cout << "------- square " << it->second.id() << std::endl;
//        std::cout << "neighbors ";
//        for (id_it=neighborIDs.begin(); id_it!=neighborIDs.end(); ++id_it) {
//            nit = _squares.find(*id_it);
//            neighbor_volume_avg += (nit->second.volume())/4.0;
//            std::cout << *id_it << ", ";
//        }
//        std::cout << std::endl;
//        // With fixed area for each square, only ever change the height
//        std::cout << "pre-smoothing volume: " << it->second.volume() << std::endl;
//        std::cout << "smoothed volume: " << neighbor_volume_avg << std::endl;
//        std::cout << "area: " << it->second.area() << std::endl;
//        std::cout << "pre-smoothing height: " << it->second.height() << std::endl;
//        std::cout << "smoothed height: " << neighbor_volume_avg/it->second.area() << std::endl;
//        // Store the result from this averaging, do not change actual height until
//        // we've computed the updated_height for every square.
//        it->second.set_updated_height(neighbor_volume_avg/it->second.area());
//    }
//    
//    // Loop again over squares to set new smoothed height
//    for (it=_squares.begin(); it!=_squares.end(); ++it) {
//        it->second.set_height(it->second.updated_height());
//    }
//    
//}


// Move the water from a Square given its current velocity and acceleration.
// Partition the volume and momentum into the neighboring Squares.
void tsunamisquares::World::moveSquares(const double dt) {
    std::map<UIndex, Square>::iterator sit;
    bool debug = true;
    
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
        SquareIDSet neighbors;
        SquareIDSet::const_iterator nit;
        
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
        
            // Find the 4 nearest squares to the new position
            neighbors = getNearestIDs(new_pos);
            
            // Init these for renormalizing the fractions
            double fraction_sum = 0.0;
            std::map<UIndex, double> originalFractions;
            std::map<UIndex, double> renormFractions;
            std::map<UIndex, double>::iterator frac_it;
            SquareIDSet to_erase;
            
            // Iterate through neighbors once to compute the fractional area overlap.
            for (nit=neighbors.begin(); nit!=neighbors.end(); ++nit) {
                // This iterator will give us the neighbor square 
                std::map<UIndex, Square>::iterator neighbor_it = _squares.find(*nit);
                double dx = fabs(new_pos[0] - squareCenter(neighbor_it->first)[0]);
                double dy = fabs(new_pos[1] - squareCenter(neighbor_it->first)[1]);
                double L = sit->second.length();
                double this_fraction = (1-dx/L)*(1-dy/L);
                
                if (debug) {
                    std::cout << "--neighbor " << *nit << std::endl;
                    std::cout << "L: " << L << std::endl;
                    std::cout << "dx: " << dx << std::endl;
                    std::cout << "dx/L: " << dx/L << std::endl;
                    std::cout << "dy: " << dy << std::endl;
                    std::cout << "dy/L: " << dy/L << std::endl;
                }
                
                //assertThrow(this_fraction < 1, "Area fraction must be less than 1.");
                if (this_fraction > 0) {
                    fraction_sum += this_fraction;
                    originalFractions.insert(std::make_pair(*nit, this_fraction));
                } else {
                    // if the fraction is less than 0, then it is not a valid neighbor
                    to_erase.insert(*nit);
                }
            }
            
            // Remove invalid neighbors from the neighbors set
            for (nit=to_erase.begin(); nit!=to_erase.end(); ++nit) {
                neighbors.erase(*nit);
            }
            
            //std::cout << "summed (over " << originalFractions.size() << ") Volume fraction: " << fraction_sum << std::endl;
            
            // Then normalize these fractions to enforce conservation.
            for (frac_it=originalFractions.begin(); frac_it!=originalFractions.end(); ++frac_it) {
//                std::cout << "---Neighbor : " << frac_it->first << std::endl;
//                std::cout << "original fraction : " << frac_it->second << std::endl; 
//                std::cout << "renormed fraction : " << (frac_it->second)/fraction_sum << std::endl; 
                //assertThrow((frac_it->second)/fraction_sum < 1, "Area fraction must be less than 1.");
                renormFractions.insert(std::make_pair(frac_it->first, (frac_it->second)/fraction_sum));
            }
            
            // Check that the normalized fractions sum exactly to 1
            double renormSum = 0.0;
            for (frac_it=renormFractions.begin(); frac_it!=renormFractions.end(); ++frac_it) {
                renormSum += frac_it->second;
            }
//            std::cout.precision(17);
//            std::cout << "renormed sum Volume fraction: " << std::fixed << renormSum << std::endl;
//            assertThrow(renormSum == 1.0, "Renormed sum does not equal 1.0");
            
            // Compute height and momentum imparted to neighbors
            for (nit=neighbors.begin(); nit!=neighbors.end(); ++nit) {
                // This iterator will give us the neighbor square 
                std::map<UIndex, Square>::iterator neighbor_it = _squares.find(*nit);
                // This iterates through the renormalized fractions
                frac_it = renormFractions.find(*nit);
                double areaFraction = frac_it->second;
                
                // Update the amount of water in the neighboring square (conserve volume)
                double dV = squareVolume(sit->first)*areaFraction;
                double H = neighbor_it->second.updated_height();
                double A_n = neighbor_it->second.area();
                neighbor_it->second.set_updated_height(H + dV/A_n);
                
                // Conserve momentum, update the velocity accordingly (at the end)
                Vec<2> dM = new_velo*areaFraction*squareMass(sit->first);
                Vec<2> M  = neighbor_it->second.updated_momentum();
                neighbor_it->second.set_updated_momentum(M+dM);
                
                if (debug) {
                    std::cout << "--- Neighbor : " << neighbor_it->second.id() << std::endl;
                    std::cout << "dV: " << dV << std::endl;
                    std::cout << "dM: " << dM << std::endl;
                    std::cout << "Area fraction: " << areaFraction << std::endl;
                }
            }
            
        } else {
            // For those squares that don't move, don't change anything.
            sit->second.set_updated_height(sit->second.height());
            sit->second.set_updated_momentum(squareMomentum(sit->first));
        }
    }
    
    // Loop again over squares to set new velocity and height from accumulated height and momentum
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        sit->second.set_height(sit->second.updated_height());
        sit->second.set_velocity(sit->second.updated_momentum()/squareMass(sit->first));
    }
    
}


tsunamisquares::Vec<2> tsunamisquares::World::getGradient(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator square_it = _squares.find(square_id);
    Vec<2> gradient, center_left, center_right, center_top, center_bottom;
    bool debug = false;
    
    // Initialize the 4 points that will be used to approximate the slopes d/dx and d/dy
    // for this square. These are the centers of the neighbor squares.
    Vec<2> center = squareCenter(square_id);
    double L = square_it->second.length();
    center_left   = Vec<2>(center[0]-L, center[1]);
    center_right  = Vec<2>(center[0]+L, center[1]);
    center_top    = Vec<2>(center[0], center[1]+L);
    center_bottom = Vec<2>(center[0], center[1]-L);
    
    // Find the squareID for these locations
    UIndex leftID = whichSquare(center_left);
    UIndex rightID = whichSquare(center_right);
    UIndex topID = whichSquare(center_top);
    UIndex bottomID = whichSquare(center_bottom);
    
    // TODO: Better boundary conditions. For now, just set no acceleration along boundary
    if (leftID == square_id || rightID == square_id || topID == square_id || bottomID == square_id) {
        gradient = Vec<2>(0.0,0.0);
    } else {
        double z_left = squareLevel(leftID);
        double z_right = squareLevel(rightID);
        double z_top = squareLevel(topID);
        double z_bottom = squareLevel(bottomID);

        Vec<2> center_L = squareCenter(leftID);
        Vec<2> center_R = squareCenter(rightID);
        Vec<2> center_T = squareCenter(topID);
        Vec<2> center_B = squareCenter(bottomID);
        
        gradient[0] = (z_right-z_left)/( center_L.dist(center_R) );
        gradient[1] = (z_top-z_bottom)/( center_T.dist(center_B) );
        
        if (debug) {
            std::cout << "square  " << square_id << std::endl;
            std::cout << "d/dx " << gradient[0] << std::endl; 
            std::cout << "d/dy " << gradient[1] << std::endl;
        }
    }
    
    return gradient;
}

void tsunamisquares::World::updateAcceleration(const UIndex &square_id) {
    std::map<UIndex, Square>::iterator square_it = _squares.find(square_id);
    Vec<2> grav_accel, friction_accel, gradient;
    double G = 9.80665; //mean gravitational acceleration at Earth's surface [NIST]
    
    // gravitational acceleration due to the slope of the water surface
    gradient = getGradient(square_id);
    grav_accel = gradient*G*(-1.0);
    
    // frictional acceleration from fluid particle interaction
    friction_accel = square_it->second.velocity()*(square_it->second.velocity().mag())*(square_it->second.friction())/(-1.0*(square_it->second.height()));
    
    // Set the acceleration
    square_it->second.set_accel(grav_accel + friction_accel);
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
// Get the square_id for each of the 4 closest squares to some location = (x,y)
tsunamisquares::SquareIDSet tsunamisquares::World::getNearestIDs(const Vec<2> &location) const {
    std::map<double, UIndex>                  square_dists;
    std::map<double, UIndex>::const_iterator  it;
    std::map<UIndex, Square>::const_iterator  sit;
    SquareIDSet                               neighbors;

    // Compute distance from "location" to the center of each square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        double square_dist = squareCenter(sit->first).dist(location);
        square_dists.insert(std::make_pair(square_dist, sit->second.id()));
    }
    
    // Grab the closest 4 squares and return their IDs
    for (it=square_dists.begin(); it!=square_dists.end(); ++it) {
        neighbors.insert(it->second);
        if (neighbors.size() == 4) break;
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

// Get the square_id for each of the 4 closest squares to some square square_id
tsunamisquares::SquareIDSet tsunamisquares::World::getNeighborIDs(const UIndex &square_id) const {
    std::map<double, UIndex>                  square_dists;
    std::map<double, UIndex>::const_iterator  it;
    std::map<UIndex, Square>::const_iterator  sit;
    SquareIDSet                               neighbors;
    std::map<UIndex, Square>::const_iterator  this_sit = _squares.find(square_id);

    // Compute distance from center of the input square to the center of each other square.
    // Since we use a map, the distances will be ordered since they are the keys
    for (sit=_squares.begin(); sit!=_squares.end(); ++sit) {
        // Skip the square whose neighbors we're trying to find 
        if (sit->second.id() != square_id) {
            double square_dist = squareCenter(sit->first).dist(squareCenter(this_sit->first));
            square_dists.insert(std::make_pair(square_dist, sit->second.id()));
        }
    }
    
    // Grab the closest 4 squares and return their IDs
    for (it=square_dists.begin(); it!=square_dists.end(); ++it) {
        neighbors.insert(it->second);
        if (neighbors.size() == 4) break;
    }
    
    return neighbors;
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

double tsunamisquares::World::squareMass(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // area x height x density
    return (sit->second.area())*(sit->second.height())*(sit->second.density());
}

double tsunamisquares::World::squareVolume(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // area x height
    return (sit->second.area())*(sit->second.height());
}

tsunamisquares::Vec<2> tsunamisquares::World::squareMomentum(const UIndex &square_id) const {
    std::map<UIndex, Square>::const_iterator  sit = _squares.find(square_id);
    // (x,y) of square center
    return (sit->second.velocity())*squareMass(square_id);
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
        std::cout << "mass: " << squareMass(square_id) << std::endl;
        std::cout << "velocity: " << this_square.velocity() << std::endl; 
        std::cout << "accel: " << this_square.accel() << std::endl;    
        std::cout << "momentum: " << squareMomentum(square_id) << std::endl; 
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

    out_stream << time << "\t";

    //
    for (i=0; i<2; ++i) {
        out_stream << squareLatLon(square_id)[i] << "\t\t";
    }

    out_stream << squareLevel(square_id) << "\t\t";

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
        double A, B;
        int row = (int)(i/num_lons);
        int col = (int)(i%num_lons);
        if (col==num_lons-1 && row!=num_lats-1) {
            // right edge, not bottom row
            A = (_vertices[i].xy() - _vertices[i-1].xy()).mag();
            B = (_vertices[i].xy() - _vertices[i+num_lons].xy()).mag();
        } else if (row==num_lats-1 && col!=num_lons-1) {
            // bottom row not right edge
            A = (_vertices[i].xy() - _vertices[i-num_lons].xy()).mag();
            B = (_vertices[i].xy() - _vertices[i+1].xy()).mag();
        } else if (row==num_lats-1 && col==num_lons-1) {
            // bottom right corner
            A = (_vertices[i].xy() - _vertices[i-num_lons].xy()).mag();
            B = (_vertices[i].xy() - _vertices[i-1].xy()).mag();
        } else {
            A = (_vertices[i].xy() - _vertices[i+1].xy()).mag();
            B = (_vertices[i].xy() - _vertices[i+num_lons].xy()).mag();
        }
        new_square.set_area(A*B);

//        // TEMP FIX TO FORCE REGULAR GRID
//        if (i==0) {
//            double A = (_vertices[i].xy() - _vertices[i+1].xy()).mag();
//            double B = (_vertices[i].xy() - _vertices[i+num_lons].xy()).mag();
//            new_square.set_area(A*B);
//        } else {
//            new_square.set_area(_squares[0].area());
//        }

        
        _squares.insert(std::make_pair(new_square.id(), new_square));
    }
    

    return 0;
}


int tsunamisquares::World::write_file_kml(const std::string &file_name) {
    std::ofstream                             out_file;
    std::map<UIndex, Square>::const_iterator  sit;
    LatLonDepth                               min_bound, max_bound, center;
    Vec<3>                                    min_xyz, max_xyz;
    double                                    dx, dy, range, L;
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
        L           = sit->second.length();
        // Locate the corners in XYZ, then convert to LLD
        v3      = Vec<3>(centerXY[0]-L/2.0, centerXY[1]+L/2, 0.0); // top left
        lld[0]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]-L/2.0, centerXY[1]-L/2, 0.0); // bottom left
        lld[1]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+L/2.0, centerXY[1]-L/2, 0.0); // bottom right
        lld[2]  = c.convert2LatLon(v3);
        v3      = Vec<3>(centerXY[0]+L/2.0, centerXY[1]+L/2, 0.0); // top left
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


