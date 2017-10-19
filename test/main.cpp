#include "../ldpc/ldpc_encode.hpp"
#include "../ldpc/ldpc_decode.hpp"
#include "../ldpc/ldpc_matrix.hpp"

#include "../ldpc/ldpc_build_matrix_v02.hpp"    //is necessary to build a matrix
#include "../ldpc/ldpc_utils.hpp"    //print utility


static void test_ldpc1()
{
        ldpc::encode_t<ldpc::c_num_h_columns, ldpc::c_message_bits>
        enc(ldpc::c_g_matrix, ldpc::c_g_matrix_column_order);
        ldpc::decode_soft_t<ldpc::c_num_h_columns, ldpc::c_message_bits, ldpc::c_sparsed_max_num_columns>
        dec(ldpc::c_h_matrix, ldpc::c_g_matrix_column_order);
        
        uint8_t const msg[enc.c_message_bytes] = {'B','l','a','c','k','C', 'a', 't'};
        uint8_t restored_msg[enc.c_message_bytes];
        uint8_t restored_msg_soft[enc.c_message_bytes];
        uint8_t codeword[enc.c_codeword_bytes];
        
        enc.encode(codeword, msg);
        
        bool codeword_bits_orig[ldpc::c_num_h_columns];
        bool codeword_bits[ldpc::c_num_h_columns];
        ldpc::bit_ops::convert_from_bytes(codeword_bits_orig, codeword);
        memcpy(codeword_bits, codeword_bits_orig, sizeof(codeword_bits));
        
        float factors[ldpc::c_num_h_columns];
        for (int i = 0; i < ldpc::c_num_h_columns; ++i) {
                if (codeword_bits[i])
                        factors[i] = 1.f;
                else
                        factors[i] = 0.f;
        }
        
        //now add the distortion
        codeword_bits[1] = !codeword_bits[1];
        codeword_bits[77] = !codeword_bits[77];
        codeword_bits[30] = !codeword_bits[30];
        //codeword_bits[62] = !codeword_bits[62];               //this matrix doesn't like bit 62... really!
        codeword_bits[64] = !codeword_bits[64];
        ldpc::bit_ops::convert_from_bit_array(codeword, codeword_bits);
        
        bool const ok = dec.decode_bsc_hard(restored_msg, codeword, 50);
        
        //add floating-point distortion...
        if (codeword_bits_orig[1]) factors[1] = 0.3f; else factors[1] = 0.7f;
        if (codeword_bits_orig[77]) factors[77] = 0.2f; else factors[77] = 0.6f;
        if (codeword_bits_orig[83]) factors[83] = 0.2f; else factors[83] = 0.6f;
        if (codeword_bits_orig[119]) factors[119] = 0.2f; else factors[119] = 0.6f;
        bool const soft_ok = dec.decode_bsc_soft(restored_msg_soft, factors);
        
        
        printf("test_ldpc1: hard: %d  soft: %d\n", ok, soft_ok);
}

static void build_matrix()
{
        bool mx[ldpc::v02::c_m][ldpc::v02::c_n];
        ldpc::v02::gen_H_matrix(mx);
        ldpc::print_matrix(mx);
        
        auto print_loop_stat = [&mx](ldpc::v02::loops_t const & ls) -> void
        {
                for(int i = 0; i < sizeof(ls.num)/sizeof(uint16_t); ++i)
                        printf("%d\t", ls.num[i]);
                printf("\n");
        };
        
        
        ldpc::v02::remove_loops(mx, print_loop_stat);
        ldpc::print_matrix(mx);
        ldpc::print_matrix_c_sparsed(mx);
        ldpc::print_matrix_for_make_pchk(mx);
        
        //ldpc::print_matrix(e);
}

int main() {
        test_ldpc1();
        return 0;
}
