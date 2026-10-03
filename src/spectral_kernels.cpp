#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Batched short-time Fourier transform power spectrogram.
//
// Replaces the per-window R loop in spectrogram() with a compiled radix-2 FFT
// over all analysis windows, parallelized across windows with OpenMP. The
// output matches the R reference exactly (same one-sided PSD normalization):
//   psd = |FFT(segment .* window)|^2 / (sr * window_power)
//   psd[bins 1..nf-2] *= 2                  (one-sided, non-DC/non-Nyquist)
// where nf = nfft/2 + 1 and nfft = length(window) must be a power of two
// (the R wrapper falls back to a batched mvfft path otherwise).

static inline bool is_pow2(int n) { return n > 0 && (n & (n - 1)) == 0; }

// [[Rcpp::export]]
NumericMatrix cpp_stft_power(NumericVector signal, NumericVector window,
                             int step, double sr, double window_power) {
  const int nfft = static_cast<int>(window.size());
  if (!is_pow2(nfft)) stop("cpp_stft_power requires a power-of-two window.");
  if (step < 1) stop("cpp_stft_power requires step >= 1.");
  const int nf = nfft / 2 + 1;
  const R_xlen_t n = signal.size();
  // Signal shorter than one window yields no frames. Guard before the division
  // (C++ integer division truncates toward zero, not floor, so a negative
  // (n - nfft) must not be allowed to produce a phantom window).
  if (n < nfft) return NumericMatrix(nf, 0);
  const int n_windows = static_cast<int>((n - nfft) / step + 1);
  if (n_windows <= 0) return NumericMatrix(nf, 0);

  // Bit-reversal permutation table.
  std::vector<int> rev(nfft, 0);
  int logn = 0;
  while ((1 << logn) < nfft) ++logn;
  for (int i = 0; i < nfft; ++i) {
    int x = i, r = 0;
    for (int b = 0; b < logn; ++b) { r = (r << 1) | (x & 1); x >>= 1; }
    rev[i] = r;
  }
  // Twiddle table W[k] = exp(-2*pi*i*k/nfft), k = 0 .. nfft/2 - 1.
  std::vector<double> Wr(nfft / 2), Wi(nfft / 2);
  const double twopi = 6.28318530717958647692;
  for (int k = 0; k < nfft / 2; ++k) {
    Wr[k] = std::cos(-twopi * k / nfft);
    Wi[k] = std::sin(-twopi * k / nfft);
  }

  NumericMatrix out(nf, n_windows);
  const double* sig = REAL(signal);
  const double* win = REAL(window);
  const double norm = sr * window_power;
  double* op = REAL(out);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector<double> re(nfft), im(nfft);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int w = 0; w < n_windows; ++w) {
      const R_xlen_t start = static_cast<R_xlen_t>(w) * step;
      // Load windowed segment directly into bit-reversed positions.
      for (int i = 0; i < nfft; ++i) {
        re[rev[i]] = sig[start + i] * win[i];
        im[rev[i]] = 0.0;
      }
      // Iterative radix-2 Cooley-Tukey (in place).
      for (int len = 2; len <= nfft; len <<= 1) {
        const int half = len >> 1;
        const int tstep = nfft / len;
        for (int base = 0; base < nfft; base += len) {
          for (int j = 0; j < half; ++j) {
            const int t = j * tstep;
            const double wr = Wr[t], wi = Wi[t];
            const int a = base + j, b = base + j + half;
            const double vr = re[b] * wr - im[b] * wi;
            const double vi = re[b] * wi + im[b] * wr;
            const double ur = re[a], ui = im[a];
            re[a] = ur + vr; im[a] = ui + vi;
            re[b] = ur - vr; im[b] = ui - vi;
          }
        }
      }
      // One-sided PSD.
      double* col = op + static_cast<R_xlen_t>(w) * nf;
      for (int k = 0; k < nf; ++k) {
        double p = (re[k] * re[k] + im[k] * im[k]) / norm;
        if (k > 0 && k < nf - 1) p *= 2.0;
        col[k] = p;
      }
    }
  }
  return out;
}
