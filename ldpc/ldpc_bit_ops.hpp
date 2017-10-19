#pragma once
#include <stdint.h>
#include <memory.h>
#include <cmath>
#include <algorithm>

namespace ldpc { namespace bit_ops{

//mad - multiply and add
template<int c_n>
bool mad_rows(bool const (&a)[c_n], bool const (&b)[c_n]);
       
template<int c_m, int c_n>
bool mad_row_col(bool const (&r)[c_n], bool const (&c)[c_m][c_n], int const column_idx);
        
template<int c_max_num_col, int c_n>
bool mad_rows_sparse(int const (&a)[c_max_num_col], bool const (&b)[c_n]);
       
template<int c_n>
void column_reorder_enc(bool (&a)[c_n], int const (&column_idx)[c_n]);

template<int c_n>
void column_reorder_dec(bool (&a)[c_n], int const (&column_idx)[c_n]);

//returns a/b aligned to up
constexpr int align_division_up(int const a, int const b);
        
//converts any data to boolean array;
//if data size in bits is longer than c_n, the least significant bits are moved to r
template<int c_n>
void convert_from_bytes(bool (&r)[c_n], uint8_t const (&data)[align_division_up(c_n,8)]);

//converts boolean array to data back. if there is bit-space in r, most significant beats are filled with zero
template<int c_n>
void convert_from_bit_array(uint8_t (&r)[align_division_up(c_n,8)], bool const (&)[c_n]);

//llr - log likelyhood 
template <int c_n>
void rasterize_llr(bool (&codeword)[c_n], float const (&probs)[c_n]);

//convert probabilties of bit to be 1 to bits by 0.5 threshold
template <int c_n>
void rasterize_probs(bool (&codeword)[c_n], float const (&probs)[c_n]);

        
}}


#include "ldpc_bit_ops.impl.hpp"
