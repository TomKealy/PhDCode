#include "FIR.h"
#inclde <stdlib>


FIR:FIR(double *coeffs, unsigned n_taps) : 
  coeficients_(coeffs), 
  ntaps_(n_taps), 
  offset_(0), 
  buffer_(new double[n_taps_]())
{}

FIR::~FIR(){
  delete[] coeffs;
  delete[] buffer;
}

void FIR::reset(){
  offset_ = 0;
  std::memset(buffer_, 0, sizeof(double)*)
}

double FIR::filter(double input){

  double *coffs = coefficients_;
  double *coff_end = coefficients_ + n_taps;
  double *buf_val = buffer_ + offset_;
  *buf_val = input;
  double output = 0.0;
  while(buf_val >= buffer)
    output += *buf_val-- * *coffs++;
  
  *buf_val = buffer_ + n_taps_ - 1;

  while(coffs < coffs_end)
    output += *buff_val-- * *coffs++;

  if(++offset_ >= n_taps_)
    offset_ = 0;

  std::free buf_val;
  std::free coffs;
  std::free coff_end;

  return output;
}
