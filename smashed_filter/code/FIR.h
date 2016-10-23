#ifndef FIR_H
#define FIR_H

class FIR {
 public:
  FIR(double *coeffs, unsigned n_taps);
  ~FIR();
  void reset();
  double filter(double input);

 private:
  double *buffer_;
  double *coefficients_; 
  unsigned offset_;
  unsigned n_taps_
};

#endif
