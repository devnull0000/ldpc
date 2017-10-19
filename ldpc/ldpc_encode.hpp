#pragma once
//Yuriy Shevyrov's implementation of LDPC codes; Sep 11 2017
//see http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/progs.html [site of Radford M. Neal + his software, great thanks to him!]
//and http://sigpromu.org/sarah/SJohnsonLDPCintro.pdf  [Introducing Low-Density Parity-Check Codes; Sarah J. Johnson]

#include "ldpc_bit_ops.hpp"
#include <vector>
#include <memory.h>

namespace ldpc {
        
//c_n - total number of bits in codeword
//c_k - length of an original message to encode
template<int c_n, int c_k>
struct encode_t
{
        constexpr static int c_check_bits = c_n - c_k;
        constexpr static int c_codeword_bytes = bit_ops::align_division_up(c_n,8);
        constexpr static int c_message_bytes = bit_ops::align_division_up(c_k,8);
        
        //from where to take the matrix ? read ldpc_matrix.hpp or generate it yourself!
        //matrix is assumed to be in constant memory (and be valid during encoder_t work)
        encode_t(bool const (&g_matrix)[c_check_bits][c_k], int const (&column_reorder_idx)[c_n]) : m_g(g_matrix), m_column_idx(column_reorder_idx){}
        
        void                 encode(bool (&cw)[c_n], bool const (&msg)[c_k]);
        void                 encode(uint8_t (&r)[c_codeword_bytes], uint8_t const (&msg)[c_message_bytes]);
        
        //it encodes all the data stream, just a for-loop wrapper of the encode function above
        //it assumes that both input and output bits (c_n, c_k) can be divided by 8, else it will work but will not compress them in a single bit-stream keeping gaps between packets
        std::vector<uint8_t> encode(std::vector<uint8_t> const & data);
        
private:
        bool const (&m_g)[c_check_bits][c_k];
        int  const (&m_column_idx)[c_n];
};
        
}

#include "ldpc_encode.impl.hpp"

