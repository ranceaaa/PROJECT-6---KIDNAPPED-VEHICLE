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
#define small 0.000001
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 80;  // TODO: Set the number of particles
  std::default_random_engine gen;
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  std::normal_distribution<double> dist_x(x,std_x);
  std::normal_distribution<double> dist_y(y,std_y);
  std::normal_distribution<double> dist_theta(theta,std_theta);
  
  for(int i = 0;i<num_particles;i++){
  
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);

  }
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
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  std::default_random_engine gen;

  std::normal_distribution<double> dist_x(0,std_x);
  std::normal_distribution<double> dist_y(0,std_y);
  std::normal_distribution<double> dist_theta(0,std_theta);
 
  for(int i = 0; i<num_particles;i++)
  {
    double theta = particles[i].theta;
    if(fabs(yaw_rate) < small){
	particles[i].x += velocity*delta_t*cos(theta);
        particles[i].y += velocity*delta_t*sin(theta);
    }else{
	particles[i].x += (velocity/yaw_rate)*(sin(theta+yaw_rate*delta_t) - sin(theta));
        particles[i].y += (velocity/yaw_rate)*(cos(theta) - cos(theta+yaw_rate*delta_t));
        particles[i].theta += yaw_rate * delta_t;
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
    
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
  int num_predictions = predicted.size();
  int num_observations = observations.size();
  
  for(int i = 0; i < num_observations; i++)
  {
    double min_dist = 99999999.99;
    int id = -99;
    for(int j = 0; j < num_predictions; j++)
    {
      double x_distance = observations[i].x - predicted[j].x;
      double y_distance = observations[i].y - predicted[j].y;
      double distance = x_distance*x_distance + y_distance*y_distance;

      if(distance < min_dist)
      {
	min_dist = distance;
        id = predicted[j].id;
      } 

    }
    observations[i].id = id;
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
  double std_range = std_landmark[0];
  double std_bearing = std_landmark[1];
  for(int i = 0; i<num_particles;i++)
  {
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    double sensor_range_la2 = sensor_range*sensor_range;
    vector<LandmarkObs> rangeLandmarks;
    for(unsigned int j = 0; j<map_landmarks.landmark_list.size(); j++)
    {
      double landmark_x = map_landmarks.landmark_list[j].x_f;
      double landmark_y = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dif_x = x - landmark_x;
      double dif_y = y - landmark_y;
      if(dif_x*dif_x+dif_y*dif_y <= sensor_range_la2)
      {
	rangeLandmarks.push_back(LandmarkObs{id, landmark_x, landmark_y});
      }
    }
    vector<LandmarkObs> map_obs;
    
    for(unsigned int j = 0; j<observations.size();j++)
    {
      double x_map = x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
      double y_map = y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
      map_obs.push_back(LandmarkObs{ observations[j].id, x_map, y_map});
    }
    dataAssociation(rangeLandmarks,map_obs);
    int num_range = rangeLandmarks.size();
    double weight = 1.0;
    for(unsigned int j = 0; j<map_obs.size();j++)
    {
      bool cont = false;
      int k = 0;
      double new_x = 0.0;
      double new_y = 0.0;
      while(cont==false && k < num_range)
      { 
        if(rangeLandmarks[k].id == map_obs[j].id)
        {
	  cont = true;
          new_x = rangeLandmarks[k].x;
          new_y = rangeLandmarks[k].y;
        }
        k += 1;
      }
      double mu_x = map_obs[j].x - new_x;
      double mu_y = map_obs[j].y - new_y;
      double gauss_norm;
      gauss_norm = 1 / (2 * M_PI * std_range * std_bearing);

  // calculate exponent
      double exponent;
      exponent = (pow(mu_x, 2) / (2 * pow(std_range, 2)))+ (pow(mu_y, 2) / (2 * pow(std_bearing, 2)));
    
  // calculate weight using normalization terms and exponent
      double weigh= gauss_norm * exp(-exponent);
      if(weigh == 0)
      {
        weigh = 0.0001;
      }
      weight *= weigh;
      particles[i].weight = weight;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> weights;
  for(int i = 0; i< num_particles;i++)
  { double part_weight = particles[i].weight;
    weights.push_back(part_weight);
  }
  // I was stucked for a while at this point and the question with number 313824 saved my day
  std::default_random_engine gen;
  std::discrete_distribution<> dist(weights.begin(), weights.end());
  
  std::vector<Particle> parti;
  while (parti.size() < particles.size()) {
    int id = dist(gen);
    parti.push_back(particles[id]);
  }
  
  particles = parti;
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
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
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
