/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  particles.resize(num_particles);
  weights.resize(num_particles, 1.0);
  // create normal (Gaussian) distributions for x,y,theta
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  // initialize all particles in loops
  for (int i = 0; i < num_particles; ++i) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;

  }

  // set initilization flag to true
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   // create normal (Gaussian) distributions for x,y,theta
   std::default_random_engine gen;
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);

   for (int i=0; i< num_particles; ++i) {
     if (fabs(yaw_rate) <1e-4) {
        particles[i].x += velocity * delta_t * cos(particles[i].theta);
        particles[i].x += dist_x(gen);
        particles[i].y += velocity * delta_t * sin(particles[i].theta);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
     }
     else{
        particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta+yaw_rate * delta_t) - sin(particles[i].theta));
        particles[i].x += dist_x(gen);
        particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta+yaw_rate * delta_t));
        particles[i].y += dist_y(gen);
        particles[i].theta += yaw_rate * delta_t;
        particles[i].theta += dist_theta(gen);
     }

     // normalize theta
     // particles[i].theta = std::fmod(particles[i].theta, 2.0 * M_PI);

   }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   for (unsigned int i=0; i<observations.size(); ++i){
      const auto& x = observations[i].x;
      const auto& y = observations[i].y;
      observations[i].id = -1;       // set the observation id to invalid value first
      double thres = numeric_limits<double>::max();     // initialise threshold value to infinite
      for (unsigned int j=0; j<predicted.size(); ++j) {
        double distance = dist(x,y,predicted[j].x,predicted[j].y);
        if (distance < thres) {
            thres = distance;        // update threshold
            observations[i].id = predicted[j].id;       // use observations id property to keep nearest landmark id
        }
      }
   }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

   for (int i = 0; i<num_particles; ++i) {
       // initialise vector of observation in map coordinates, size equals to observation.
       vector<LandmarkObs> observations_MAP(observations.size());
       // initialise vector of landmarks for current particle.
       vector<LandmarkObs> landmarks;     // empty
       // initialise a temp variation to preserve current particle.
       const auto& particle = particles[i];
       // transform observation from Vehicle's coordinate location to MAP's coordinate location
       vcl2MAP(particle,observations,observations_MAP);
       // get current particle's landmark within sensor_range
       MAP2landmark(particle,sensor_range,map_landmarks,landmarks);
       // call dataAssociation to match landmarks to the nearest observations_MAP point
       dataAssociation(landmarks,observations_MAP);

       vector<int> associations(observations_MAP.size());
       vector<double> sense_x(observations_MAP.size());
       vector<double> sense_y(observations_MAP.size());
       for (unsigned int i=0; i<observations_MAP.size(); ++i){
          associations[i] = observations_MAP[i].id;
          sense_x[i] = observations_MAP[i].x;
          sense_y[i] = observations_MAP[i].y;
       }
       SetAssociations(particles[i], associations, sense_x, sense_y);
       // compute each particle's weight based on Multivariate_normal_distribution
       // particles[i].weight = 1.0;
       double prob = 1.0;
       bool found = false;
       for (unsigned int i=0; i<observations_MAP.size(); ++i) {
          // cout<<"size of observations_MAP: "<<observations_MAP.size()<<endl;
          for (unsigned j=0; j<landmarks.size(); ++j){
             // cout<<"i: "<<i<<", "<<"j: "<<j<<endl;
             if (landmarks[j].id == observations_MAP[i].id){
                found = true;
                // cout<<"found landmark["<<j<<"] match with obeservation["<<i<<"]"<<endl;
                double mu_x = landmarks[j].x;
                // cout<<"mu_x: "<<mu_x<<endl;
                double mu_y = landmarks[j].y;
                // cout<<"mu_y: "<<mu_y<<endl;
                prob *= multiv_prob(std_landmark[0], std_landmark[1], observations_MAP[i].x, observations_MAP[i].y,
                   mu_x, mu_y);
                // cout<<"weight updated: "<<prob<<endl;
                // break;
             }
          }

       }
       if(!found){
          particles[i].weight = 0.0;
       }
       else{
          particles[i].weight = prob;
       }

       // cout<<"particles["<<i<<"].weight:"<<particles[i].weight<<endl;
   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   double total_weights = 0.0;
   for (int i=0; i<num_particles; ++i){
      total_weights += particles[i].weight;
   }

   // weights vector is used for initialization of discrete_distribution
   for (int i=0; i<num_particles; ++i){
      weights[i] = (particles[i].weight/total_weights);
   }

   /*
   // print weights
   for (unsigned int i=0; i<weights.size(); ++i){
      cout<<weights[i]<<" ";
   }
   */

   // create a temp vector of particles for resample, after sampling,
   // replace the origin with resampled particles vector
   // vector<Particle> p1(num_particles);
   vector<Particle> p1;

   std::default_random_engine gen;
   std::discrete_distribution<int> distribution (weights.begin(),weights.end());
   for (int i=0; i<num_particles; ++i){
      int index = num_particles;      // set index to invalid value;
      // loop until that generate index within particles range
      while(index >= num_particles){
         index = distribution(gen);
      }
      // insert the sampled particle into empty p1
      p1.push_back(particles[index]);
      // p1[i] = particles[index];
   }
   // cout<<"size of p1: "<<p1.size()<<endl;
   // assign p1 to particles
   particles.assign(p1.begin(),p1.end());
   // particles.swap(p1);

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

void ParticleFilter::vcl2MAP(const Particle& particle, const std::vector<LandmarkObs> & observations, std::vector<LandmarkObs> & observations_MAP)
{
    for (unsigned int i=0; i<observations.size(); ++i){
        observations_MAP[i].x = particle.x + (cos(particle.theta) * observations[i].x - sin(particle.theta) * observations[i].y);
        observations_MAP[i].y = particle.y + (sin(particle.theta) * observations[i].x + cos(particle.theta) * observations[i].y);
    }
}

void ParticleFilter::MAP2landmark(const Particle& particle,double sensor_range,const Map &map_landmarks,vector<LandmarkObs> & landmarks)
{
    for (unsigned int i=0; i<map_landmarks.landmark_list.size(); ++i) {
        if (dist(particle.x,particle.y,map_landmarks.landmark_list[i].x_f,map_landmarks.landmark_list[i].y_f) <sensor_range){
            LandmarkObs obs;
            obs.id = map_landmarks.landmark_list[i].id_i;
            obs.x = map_landmarks.landmark_list[i].x_f;
            obs.y = map_landmarks.landmark_list[i].y_f;
            landmarks.push_back(obs);
        }
    }
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);

  return weight;
}