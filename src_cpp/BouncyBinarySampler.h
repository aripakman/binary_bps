/* 
 * File:   HMC_BinarySampler.h
 * Author: aripakman
 *
 * Created on May 8, 2013, 10:37 PM
 */

#ifndef HMC_BINARYSAMPLER_H
#define	HMC_BINARYSAMPLER_H

#include <vector>
#include "BinaryDistribution.h"


using namespace std;
typedef pair< vector<int> ,double> ST_Pair;    

class BouncyBinarySampler {
    
public:
   
    BouncyBinarySampler(const BinaryDistribution * b ): bin_dist(b) {}     // b is passed as a pointer because an object of type BinaryDistribution 
                                                                         // cannot be created, since it contains pure virtual functions 
    
    void runSamplerGaussian(const double & max_time, vector< ST_Pair > & STs, vector<double> & log_likes, int * seed_ptr, vector<double> & Y );
    void runSamplerExponential(const double & max_time, vector< ST_Pair > & STs, vector<double> & log_likes, int * seed_ptr, vector<double> & Y );

    int get_dim() {return bin_dist->getDimension();}
    
private:
    
    
    const BinaryDistribution * bin_dist;

};

#endif	/* HMC_BINARYSAMPLER_H */

