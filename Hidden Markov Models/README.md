# Hidden Markov Models

Here I model the behavior of a robot that walks around a ring. The ring is divided into 10 sectors. At any given time point, the robot is in one of the sectors and decides with equal probability to stay in that sector or move to the next sector. The robot is equipped with a tracking device. The device is not very accurate though: If the robot is in the sector i, then the device will report that the robot is in the sectors [i − 2; i + 2] with equal probability.

## Tasks:

1. Build a hidden Markov model
2. Simulate the HMM for 100 time steps
3. Discard the hidden states from the sample obtained above. Use the remaining observations to compute the filtered and smoothed probability distributions for each of the 100 time points. Compute also the most probable path.
4. Compute the accuracy of the filtered and smoothed probability distributions, and of the most probable path.
5. Repeat the previous task with different simulated samples.
6. Does more observations result in better knowledge of where the robot is?
7. Consider any of the samples above of length 100. Compute the probabilities of the hidden states for the time step 101.
