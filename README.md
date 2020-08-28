# Kidnapped Vehicle

Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

[//]: # (Image References)

[image1]: ./images/accuracy.png "Accuracy"

---

## Prerequisites

* Term 2 Simulator with uWebSocketIO
* cmake >= 3.5
* make >= 4.1
* gcc/g++ >= 5.4

If you want to install the libraries on linux Operating System use install-linux.sh(./install-linux.sh) or for OSX use install-mac.sh(./install-mac.sh).

---

## Basic Build Instructions

1. Clone this repo.
2. Select build directory: `cd build`
3. Compile: `cmake .. && make` 
4. Run it: `./particle_filter `

---

## Results

![alt text][image1]

As can be seen from the pictures, the filter meets the required accuracy and performance.
