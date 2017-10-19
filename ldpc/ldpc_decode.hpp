#pragma once
//Yuriy Shevyrov's implementation of LDPC codes; Sep 11 2017
//see http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/progs.html [site of Radford M. Neal + his software, great thanks to him!]
//and http://sigpromu.org/sarah/SJohnsonLDPCintro.pdf  [Introducing Low-Density Parity-Check Codes; Sarah J. Johnson]

#include "ldpc_bit_ops.hpp"
#include "ldpc_tanner_graph.hpp"

//H-matrix matrix for decoding/checking result (aka syndrome)
//if syndrome iz 0 there is no error

//c_n - total length of the code
//c_k - number of message bits.
//c_max_hc - maximum number of columns in sparsed H matrix; number in a row shows where is 1. -1 indicates end of the row

namespace ldpc {

template<int c_n, int c_k, int c_max_hc>
struct decode_t
{
        constexpr static int c_check_bits = c_n - c_k;
        constexpr static int c_codeword_bytes = bit_ops::align_division_up(c_n,8);
        constexpr static int c_message_bytes = bit_ops::align_division_up(c_k,8);               //i.e. original message without parity bits
        
        decode_t(int const (&sparsed_h_matrix)[c_check_bits][c_max_hc], int const (&column_reorder_idx)[c_n]) : m_h(sparsed_h_matrix), m_column_idx(column_reorder_idx){}

        //decode corrupted data on BSC (Binary Symmetric Channel) without providing any idea of what bits were corrupted
        //returns true if the decoder thinks that the message was decoded successfully
        //the algorithm is iterative, so you should specify maximum number of iterations before giving up
        bool decode_bsc_hard(bool (&msg)[c_k], bool const (&codeword)[c_n], int const num_iterations);        
        bool decode_bsc_hard(uint8_t (&msg)[c_message_bytes], uint8_t (&codeword)[c_codeword_bytes], int const num_iterations = 50);
        
        
        bool check_codeword(bool const (&codeword)[c_n]);
        
protected:
        int const (&m_h)[c_check_bits][c_max_hc];
        int const (&m_column_idx)[c_n];
};


//********************************

template<int c_n, int c_k, int c_max_hc>
struct decode_soft_t : public decode_t<c_n, c_k, c_max_hc>
{
        constexpr static int c_check_bits = c_n - c_k;
        constexpr static int c_codeword_bytes = bit_ops::align_division_up(c_n,8);
        constexpr static int c_message_bytes = bit_ops::align_division_up(c_k,8);               //i.e. original message without parity bits

        decode_soft_t(int const (&sparsed_h_matrix)[c_check_bits][c_max_hc], int const (&column_reorder_idx)[c_n]) :
                decode_t<c_n, c_k, c_max_hc>(sparsed_h_matrix, column_reorder_idx), m_graph(sparsed_h_matrix){}

        bool decode_bsc_soft(bool (&msg)[c_k], float const (&probs)[c_n], int const num_iterations);
        bool decode_bsc_soft(uint8_t (&msg)[c_message_bytes], float const (&probs)[c_n], int const num_iterations = 50);

        bool check_codeword_and_write(bool (&msg)[c_k], float const (&probs)[c_n]);

private:
        typedef tanner_graph_t::check_node_t    check_node_t;
        typedef tanner_graph_t::bit_node_t      bit_node_t;

        void compute_out_estimates(check_node_t * const cn);    //out if look from a check node
        void compute_in_estimates(bit_node_t * const bn);
        bool check_estimates_llr(bool (&msg)[c_k]);    //and write message if everything is ok
        bool check_codeword_and_write_llr(bool (&msg)[c_k], float const (&llrs)[c_n]);


        tanner_graph_t          m_graph;
};


} //ldpc

#include "ldpc_decode.impl.hpp"
