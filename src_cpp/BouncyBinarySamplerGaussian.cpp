/* 
 * Author: aripakman 
 * Created on Oct 17, 2015
 */

#include <cstdlib>
#include <cmath>
#include <algorithm> // inner_product, max
#include <numeric>  // std::accumulate
#include <random>
#include <iostream>
#include <vector>
#include <utility>  //for the pair templated class
#include "BinaryDistribution.h"
#include "BouncyBinarySampler.h"

using namespace std;

void move_Y(vector<double> & , int , vector<double> & , double );
void sample_velocity(vector<double> & V, int d, default_random_engine & e, normal_distribution<> & normal );

typedef pair< vector<int> ,double> ST_Pair;    

void BouncyBinarySampler::runSamplerGaussian(const double & max_time, vector< ST_Pair > & STs, vector<double> & log_likes, int * seed_ptr, vector<double> & Y ){

    cout << "Gaussian Running" << endl;

    int d = bin_dist->getDimension();

    int seed = *seed_ptr;
    
    default_random_engine e(seed);
    normal_distribution<> normal(0,1);
    uniform_real_distribution<double> uniform(0.0,1.0);


    
    STs.clear();
    log_likes.clear();

    ST_Pair ST;

    vector<double> V(d);    // velocities
    vector<int> S(d);       //signs of Y    

    // sample/set initial positions
    for (int j=0; j < d; j++ ){                        
        //Y[j] = normal(e);
        S[j] = (Y[j] < 0) ? -1 : 1;
        
    }
       

    // sample initial velocity
    sample_velocity(V,  d, e, normal);

    double elapsed_time =0;
    double u, mlu, mlu_elapsed, t_min, t_grad;
    int i_min;


    int count_sat = 0;
    int count_free = 0;
    int count_payed = 0;
    int count_not_payed = 0;
    int count_refresh =0;

    double pot_change;    
    bool just_min = false;    

    double this_t = 0;
    double ll;

    while (elapsed_time < max_time){

        mlu = -log(uniform(e));
        mlu_elapsed = 0;

        int r = 0;

        while(mlu_elapsed < mlu){

            double v_y = 0;            
            for (int i=0; i< d; i++){                          
                v_y += Y[i]*V[i];     }               
            
            // find time to saturate the u budget without hitting Y[i]=0             
            if (just_min or v_y > 0){
                t_grad = - v_y + sqrt(v_y*v_y + 2*(mlu-mlu_elapsed));
                just_min = false;
            }
            else{
                t_grad = - v_y;
            }


            // find first coordinate to reach 0 
            t_min = -1;
            for (int i=0; i<d; i++){
                if (V[i]*S[i] < 0){

                        double t_i = -Y[i]/V[i];                    
                        if (t_min == -1){
                            t_min = t_i;
                            i_min = i; }  

                        else if (t_i < t_min){
                            t_min = t_i;
                            i_min = i;                                        
                        }
                }            
            }

    

            if (t_min == -1 || (t_grad != -1  &&    t_grad < t_min) ){             // either no coordinate is headed to 0 or u will saturate first

                this_t += t_grad;                                    


                move_Y(V, d, Y, t_grad);

                if (v_y > 0) {
                    mlu_elapsed = mlu;

                    double Y2 = accumulate(Y.begin(), Y.end(), 0.0, [] (double a, double b) {return a + b*b;}); 
                    v_y = 0;            
                    for (int i=0; i< d; i++){                           //update the values of Y
                        v_y += Y[i]*V[i]; }               

                    for (int i=0; i< d; i++){                
                        V[i] = V[i] -2*v_y*Y[i]/Y2; }  
                }
                else 
                {
                    just_min = true;  }

              
            }   
            else {
                this_t += t_min;                                    
                if (v_y > 0)                                    
                    mlu_elapsed += v_y*t_min + .5*t_min*t_min;


                move_Y(V, d, Y, t_min);
                Y[i_min] = 0;                                       //for numerical stablity, make sure the coordinate reaches zero
                pot_change = S[i_min]*bin_dist->getLogLikeDifference(S,i_min);


                if (pot_change <= 0 ||  pot_change < mlu-mlu_elapsed ){                // no price to pay to cross Y[i_min] = 0 or there is a price to cross Y[i_min] = 0, but we can pay it

                    ST.first = S;                // since we will change S, store its value and the time this_t spent on this value of S
                    ST.second = this_t;
                    ll = bin_dist->getLogLike(S);                        

                    STs.push_back(ST);        
                    log_likes.push_back(ll);

                    S[i_min] = -S[i_min]; 

                    elapsed_time += this_t;        
                    this_t = 0;

                    if (pot_change> 0)
                        mlu_elapsed += pot_change;
                    }
                else{                           // we cannot pay the price to cross Y[i_min] = 0, so we stop and invert the sign of V[i_min]
                    mlu_elapsed = mlu;
                    V[i_min] = -V[i_min];
                    }
                }
        }  // while(mlu_elapsed < mlu)              
    } // while (elapsed_time < max_time){
}


inline void move_Y(vector<double> & V, int d, vector<double> & Y, double t)  {
                for (int i=0; i< d; i++){                           //update the values of Y
                    Y[i] += t*V[i];
                }               
}


inline void sample_velocity(vector<double> & V, int d, default_random_engine & e, normal_distribution<> & normal ){

        for (int j=0; j < d; j++ ){                        
                V[j] = normal(e);
            }
            // normalize so V[j] is uniformly sampled from the unit sphere. 
            // Using map-reduce style.    
            double V2 = accumulate(V.begin(), V.end(), 0.0, [] (double a, double b) {return a + b*b;}); 

            V2 = sqrt(V2);
            transform(V.begin(), V.end(), V.begin(), [V2] (double v) {return v/V2;}) ;

            //cout << "V[0]: " << V[0] << endl;


}




