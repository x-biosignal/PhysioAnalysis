#include <Rcpp.h>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Fixed-length epoch extraction.
//
// Replaces the per-channel fancy-indexing loop and the NA-initialized 4D array
// allocation in epochData()'s standard path with direct contiguous copies.
// Because time is the leading (column-major) dimension, one epoch's samples for
// a fixed (channel, sample) are contiguous in both the source and destination,
// so each is a single memcpy. Parallelized across epochs with OpenMP.
//
// data     : flattened column-major array of dims (n_time, n_channels, n_samples)
// n_time, n_channels, n_samples : source dims
// centers  : 1-based event sample positions (one per valid epoch)
// pre, post: samples before / after the event; epoch_length = pre + post + 1
// Returns a flattened column-major array of dims
//   (epoch_length, n_channels, n_epochs, n_samples).
// All (center - pre) .. (center + post) windows must lie within [1, n_time];
// the R caller guarantees this by pre-filtering to in-bounds events.
//
// [[Rcpp::export]]
NumericVector cpp_epoch_fixed(NumericVector data, int n_time, int n_channels,
                              int n_samples, IntegerVector centers,
                              int pre, int post) {
  const int L = pre + post + 1;
  if (L <= 0) stop("cpp_epoch_fixed requires pre + post + 1 > 0.");
  const R_xlen_t n_epochs = centers.size();
  // Defensive bounds check: every window [center-pre, center+post] must lie
  // within [1, n_time]. The R caller pre-filters to in-bounds events, but the
  // kernel guards independently so a direct call cannot over-read.
  for (R_xlen_t e = 0; e < n_epochs; ++e) {
    const int c = centers[e];
    if (IntegerVector::is_na(c) || c - pre < 1 || c + post > n_time) {
      stop("cpp_epoch_fixed: an epoch window falls outside the signal bounds.");
    }
  }
  const R_xlen_t out_len = static_cast<R_xlen_t>(L) * n_channels * n_epochs *
    n_samples;
  NumericVector out(out_len);
  const double* src = REAL(data);
  double* dst = REAL(out);

  const R_xlen_t src_chan_stride = n_time;                       // step to next channel
  const R_xlen_t src_samp_stride = static_cast<R_xlen_t>(n_time) * n_channels;
  const R_xlen_t dst_chan_stride = L;                            // step to next channel
  const R_xlen_t dst_epoch_stride = static_cast<R_xlen_t>(L) * n_channels;
  const R_xlen_t dst_samp_stride = dst_epoch_stride * n_epochs;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (R_xlen_t e = 0; e < n_epochs; ++e) {
    const R_xlen_t start_t = static_cast<R_xlen_t>(centers[e]) - 1 - pre; // 0-based
    for (int s = 0; s < n_samples; ++s) {
      for (int ch = 0; ch < n_channels; ++ch) {
        const double* sp = src + start_t + ch * src_chan_stride +
          static_cast<R_xlen_t>(s) * src_samp_stride;
        double* dp = dst + e * dst_epoch_stride + ch * dst_chan_stride +
          static_cast<R_xlen_t>(s) * dst_samp_stride;
        std::memcpy(dp, sp, static_cast<size_t>(L) * sizeof(double));
      }
    }
  }
  return out;
}
