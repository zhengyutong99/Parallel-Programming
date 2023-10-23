#include "PPintrin.h"

// implementation of absSerial(), but it is vectorized using PP intrinsics
void absVector(float *values, float *output, int N)
{
  __pp_vec_float x;
  __pp_vec_float result;
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {

    // All ones
    maskAll = _pp_init_ones();

    // All zeros
    maskIsNegative = _pp_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _pp_vload_float(x, values + i, maskAll); // x = values[i];

    // Set mask according to predicate
    _pp_vlt_float(maskIsNegative, x, zero, maskAll); // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _pp_vsub_float(result, zero, x, maskIsNegative); //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _pp_mask_not(maskIsNegative); // } else {

    // Execute instruction ("else" clause)
    _pp_vload_float(result, values + i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _pp_vstore_float(output + i, result, maskAll);
  }
}

void clampedExpVector(float *values, int *exponents, float *output, int N)
{
  //
  // PP STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //
  int Odd_last;
  __pp_vec_float x;
  __pp_vec_int y;

  __pp_vec_float result1 = _pp_vset_float(1.f);
  __pp_vec_float result2 = _pp_vset_float(1.f);
  __pp_vec_int count;

  __pp_vec_int zero_i = _pp_vset_int(0);
  __pp_vec_int ones_i = _pp_vset_int(1);
  __pp_vec_float zero_f = _pp_vset_float(0.f);
  __pp_vec_float ones_f = _pp_vset_float(1.f);
  __pp_vec_float nines = _pp_vset_float(9.999999f);

  __pp_mask maskOnes, maskyIsZero, maskyIsNotZero, maskResult2IsLarge, maskZero;

  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    maskOnes = _pp_init_ones();
    maskyIsZero = _pp_init_ones(0);
    maskResult2IsLarge = _pp_init_ones(0);

    _pp_vload_float(x, values + i, maskOnes); // float x = values[i];
    _pp_vload_int(y, exponents + i, maskOnes); // int y = values[i];

    _pp_veq_int(maskyIsZero, y, zero_i, maskOnes); // if (y == 0) {
    
    // _pp_vsub_float(result1, ones_f, zero_f, maskyIsZero); // output[i] = 1.f;
    _pp_vstore_float(output + i, result1, maskyIsZero);

    maskyIsNotZero = _pp_mask_not(maskyIsZero); // } else {
    maskZero = maskyIsNotZero;
    _pp_vsub_float(result2, x, zero_f, maskyIsNotZero); // float result = x;
    _pp_vsub_int(count, y, ones_i, maskyIsNotZero); // int count = y - 1;
    _pp_vgt_int(maskyIsNotZero, count, zero_i, maskyIsNotZero);
    while (_pp_cntbits(maskyIsNotZero) > 0)
    {
      _pp_vmult_float(result2, result2, x, maskyIsNotZero); // result *= x;
      _pp_vsub_int(count, count, ones_i, maskyIsNotZero); // count--;
      _pp_vgt_int(maskyIsNotZero, count, zero_i, maskyIsNotZero); // while (count > 0)
    }
    
    _pp_vgt_float(maskResult2IsLarge, result2, nines, maskZero); // if (result > 9.999999f)
    _pp_vadd_float(result2, zero_f, nines, maskResult2IsLarge); // result = 9.999999f;
    _pp_vstore_float(output + i, result2, maskZero); // output[i] = result;
    
    Odd_last = i;
  }

  __pp_mask maskOdd;
  x = _pp_vset_float(0.f);
  y = _pp_vset_int(0);
  maskOdd = _pp_init_ones(N - Odd_last);
  
  maskOnes = _pp_init_ones();
  maskyIsZero = _pp_init_ones(0);
  maskResult2IsLarge = _pp_init_ones(0);

  result1 = _pp_vset_float(1.f);
  result2 = _pp_vset_float(1.f);

  _pp_vload_float(x, values + Odd_last, maskOdd); // float x = values[i];
  _pp_vload_int(y, exponents + Odd_last, maskOdd); // int y = values[i];

  _pp_veq_int(maskyIsZero, y, zero_i, maskOdd); // if (y == 0) {
  
  // _pp_vsub_float(result1, ones_f, zero_f, maskyIsZero); // output[i] = 1.f;
  _pp_vstore_float(output + Odd_last, result1, maskyIsZero);

  maskyIsNotZero = _pp_mask_not(maskyIsZero); // } else {

  _pp_vsub_float(result2, x, zero_f, maskyIsNotZero); // float result = x;
  _pp_vsub_int(count, y, ones_i, maskyIsNotZero); // int count = y - 1;
  _pp_vgt_int(maskyIsNotZero, count, zero_i, maskyIsNotZero);
  while (_pp_cntbits(maskyIsNotZero) > 0)
  {
    _pp_vmult_float(result2, result2, x, maskyIsNotZero); // result *= x;
    _pp_vsub_int(count, count, ones_i, maskyIsNotZero); // count--;
    _pp_vgt_int(maskyIsNotZero, count, zero_i, maskyIsNotZero); // while (count > 0)
  }
  
  _pp_vgt_float(maskResult2IsLarge, result2, nines, maskOdd); // if (result > 9.999999f)
  _pp_vadd_float(result2, zero_f, nines, maskResult2IsLarge); // result = 9.999999f;
  _pp_vstore_float(output + Odd_last, result2, maskOdd); // output[i] = result;
  output[N] = 0.f;
  for (int i = N; i < N + VECTOR_WIDTH; i++)
  {
    output[i] = 0.f;
  }
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N)
{

  //
  // PP STUDENTS TODO: Implement your vectorized version of arraySumSerial here
  //

  float sum;
  __pp_vec_float x;
  __pp_vec_float result = _pp_vset_float(0.f);
  __pp_mask maskAll, maskSum;

  maskAll = _pp_init_ones();
  maskSum = _pp_init_ones(1);

  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    _pp_vload_float(x, values + i, maskAll);
    _pp_vadd_float(result, result, x, maskAll);
  }


  int vec_width = VECTOR_WIDTH;
  while (vec_width > 1)
  {
    _pp_hadd_float(result, result);
    _pp_interleave_float(result, result);
    vec_width /= 2;
  }
  
  _pp_vstore_float(&sum, result, maskSum);

  return sum;
}