# Adaptive-filter
Developing adaptive filtering in GNU Octave and Matlab.

You only need to run 'applied_in_music' in GNU-Octave or Matlab (but you need to change in code what is the music file).
I am using adaptive filtering to perform system identification, "guessing" what "unknown" filter u was applied to a known input x (from mp3 file)

Currently, there are 3 branches:

* master: FIR filter identification, written for Matlab, uses NLMS method (see https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS) )

* octave: FIR filter identification, written for Octave, laso uses NLMS method (see https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS) )

* octave_IIR_filter: IIR filter identification, written for Octave. IIR identification is not a simple extension of FIR identification (the problem does not have a single solution). It implements a Full gradient IIR LMS algorithm, which results in a approximation of the system. Results are much better when you already know the system's number of poles and number of zeroes.
