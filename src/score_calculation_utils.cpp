#include "java_compat.h"
#include <cmath>
#include <stdexcept>

// Porting from earlier Java implementation to match behavior
// see: https://tinyurl.com/java-implementation

// ---------------------------------------------------------
// 1. EXACT REPLICATION OF JAVA Z-STATISTICS
// ---------------------------------------------------------

// Helper function from ZStatistics.java
double oneMinusNormalCDFInversePLT5(double p) {
  if (p < 0) {
    stop("oneMinusNormalCDFInversePLT5 called with negative p");
  } else if (p > 0.5) {
    stop("oneMinusNormalCDFInversePLT5 called with p > 0.5");
  }

  double t = std::sqrt(-2.0 * std::log(p));
  double temp = 2.515517 + 0.802853 * t + 0.010328 * t * t;
  temp = t - temp / (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);
  return temp;
}

// [[Rcpp::export]]
NumericVector get_java_zscores(NumericVector pvals) {
  int n = pvals.size();
  NumericVector zscores(n);

  for (int i = 0; i < n; ++i) {
    double p = pvals[i];
    if (p <= 0.5) {
      if (p > 0) {
        zscores[i] = oneMinusNormalCDFInversePLT5(p);
      } else {
        zscores[i] = R_PosInf;
      }
    } else if (p < 1.0) {
      zscores[i] = -oneMinusNormalCDFInversePLT5(1.0 - p);
    } else {
      zscores[i] = R_NegInf;
    }
  }

  return zscores;
}


// ---------------------------------------------------------
// 2. MONTE CARLO BACKGROUND DISTRIBUTION (exact Java replica)
//
// IMPORTANT: z_scores MUST be supplied in Java's networkNodeList order
// (i.e. adjacency.keySet() iteration order, reconstructed on the R side by
// .java_node_order()). Collections.shuffle starts from that arrangement, so a
// different starting order with the same seed yields different permutations and
// therefore different means/stds. With the correct order this reproduces the
// Java means/stds to floating-point precision.
// ---------------------------------------------------------

// [[Rcpp::export]]
List get_java_mc_calibration(NumericVector z_scores, int trials = 2000, int seed = 1234) {
  int n = z_scores.size();

  // Internal accumulators use n+1 to map directly to subnetwork sizes
  std::vector<double> samplingScoreSums(n + 1, 0.0);
  std::vector<double> samplingScoreSquareSums(n + 1, 0.0);

  // Copy z-scores; this vector is shuffled in place each trial.
  std::vector<double> z_vec = as<std::vector<double>>(z_scores);

  JavaRandom rng(seed);

  for (int trial = 0; trial < trials; ++trial) {
    java_shuffle(z_vec, rng);

    double zSum = 0.0;
    for (int i = 0; i < n; ++i) {
      zSum += z_vec[i];
      int numberOfNodesInSubnetwork = i + 1;

      // Java's calculateScoreOfSubnetwork returns 0 for size-1; sizes > 1 use
      // the raw zSum/sqrt(k) here (normalisation is off during calibration).
      if (numberOfNodesInSubnetwork > 1) {
        double score = zSum / std::sqrt((double)numberOfNodesInSubnetwork);
        samplingScoreSums[numberOfNodesInSubnetwork] += score;
        samplingScoreSquareSums[numberOfNodesInSubnetwork] += score * score;
      }
    }
  }

  // Output vectors are size n; index i represents subnetwork size (i+1),
  // matching the 0-based lookup in greedy_expand (sc_means[new_size - 1]).
  NumericVector samplingScoreMeans(n);
  NumericVector samplingScoreStds(n);

  for (int i = 0; i < n; ++i) {
    int size = i + 1; // The actual subnetwork size this index represents

    if (size == 1) {
      samplingScoreMeans[i] = 0.0;
      // Size-1 networks score 0 and are bypassed; 1.0 guards against div-by-zero.
      samplingScoreStds[i] = 1.0;
    } else {
      double mean = samplingScoreSums[size] / trials;
      samplingScoreMeans[i] = mean;

      // var = E[x^2] - E[x]^2, matching the Java accumulation, with Java's
      // identical 1e-7 stabiliser inside the sqrt.
      double var = (samplingScoreSquareSums[size] / trials) - (mean * mean);
      samplingScoreStds[i] = std::sqrt(var + 0.0000001);
    }
  }

  return List::create(
    Named("means") = samplingScoreMeans,
    Named("stds")  = samplingScoreStds
  );
}
