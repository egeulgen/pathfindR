#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>

using namespace Rcpp;

// ---------------------------------------------------------
// 1. EXACT REPLICATION OF JAVA Z-STATISTICS
// Porting from earlier Java implementation to match behaviour
// see: https://tinyurl.com/java-implementation
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
// 2. MONTE CARLO BACKGROUND DISTRIBUTION
// ---------------------------------------------------------

// [[Rcpp::export]]
List get_java_mc_calibration(NumericVector z_scores, int trials = 2000, int seed = 42) {
  int n = z_scores.size();

  // Internal accumulators use n+1 to map directly to subnetwork sizes
  std::vector<double> samplingScoreSums(n + 1, 0.0);
  std::vector<double> samplingScoreSquareSums(n + 1, 0.0);

  std::vector<double> z_vec = as<std::vector<double>>(z_scores);
  std::mt19937 gen(seed);

  for (int trial = 0; trial < trials; ++trial) {
    std::shuffle(z_vec.begin(), z_vec.end(), gen);

    double zSum = 0.0;
    for (int i = 0; i < n; ++i) {
      zSum += z_vec[i];
      int numberOfNodesInSubnetwork = i + 1;

      if (numberOfNodesInSubnetwork > 1) {
        double score = zSum / std::sqrt((double)numberOfNodesInSubnetwork);
        samplingScoreSums[numberOfNodesInSubnetwork] += score;
        samplingScoreSquareSums[numberOfNodesInSubnetwork] += score * score;
      }
    }
  }

  // Output vectors are size n to match C++ greedy_expand 0-based lookup
  NumericVector samplingScoreMeans(n);
  NumericVector samplingScoreStds(n);

  for (int i = 0; i < n; ++i) {
    int size = i + 1; // The actual subnetwork size this index represents

    if (size == 1) {
      samplingScoreMeans[i] = 0.0;
      // Failsafe: Set std to 1.0 instead of 0.0 to strictly prevent div-by-zero
      // even though size 1 networks are bypassed in scoring.
      samplingScoreStds[i] = 1.0;
    } else {
      double mean = samplingScoreSums[size] / trials;
      samplingScoreMeans[i] = mean;

      double var = (samplingScoreSquareSums[size] / trials) - (mean * mean);
      samplingScoreStds[i] = std::sqrt(var + 0.0000001);
    }
  }

  return List::create(
    Named("means") = samplingScoreMeans,
    Named("stds")  = samplingScoreStds
  );
}
