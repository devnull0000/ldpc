#pragma once
//matrix builder as described in https://www.researchgate.net/publication/26512610_Design_LDPC_Codes_without_Cycles_of_Length_4_and_6
//unfortunately it generates only big matrices.... and I need much smaller
//so there is version 2 nearby (_v02.hpp)
namespace ldpc {
namespace v01 {

//variable names are used as in the authors' work above [so they differ from the rest of the ldpc library]
constexpr int c_j = 3;                          //number of 1s in each column
constexpr int c_k = 5;                          //number of 1s in each row
constexpr int c_n = c_k * c_k * c_k;            //number of columns; aka length of codeword (in encoded form)
constexpr int c_m = c_n * c_j / c_k;            //number of rows in the matrix & number of check nodes/bits in the message
constexpr int c_bit_nodes = c_n - c_m;          //number of message bits (i.e useful message)

constexpr int c_v = c_k;
constexpr int c_vv = c_v * c_v;
constexpr int c_vvv = c_vv * c_v;

constexpr inline void get_Dk_matrix(bool (& dk)[c_v][c_vv], int const k) {
        for( int i = 0; i < c_v; ++i )
                for( int j = 0; j < c_vv; ++j )
                        dk[i][j] = false;
        for( int i = 0; i < c_v; ++i )
                dk[i][k] = true;
}

constexpr inline void get_D_matrix(bool (& d)[c_vvv][c_vv]) {
        for( int i = 0; i < c_vv; ++i ) {
                int const wi = i * c_v;
                bool dk[c_v][c_vv] = {false};
                get_Dk_matrix(dk, i);
                for( int si = 0; si < c_v; ++si )
                        for( int sj = 0; sj < c_vv; ++sj )
                                d[wi + si][sj] = dk[si][sj];
        }
}

constexpr inline void get_Ek_matrix(bool (& ek)[c_v][c_vv], int const k) {
        for( int i = 0; i < c_v; ++i )
                for( int j = 0; j < c_vv; ++j )
                        ek[i][j] = false;
        for( int i = 0; i < c_v; ++i )
                ek[i][(i + k) % c_vv] = true;
}

constexpr inline void get_E_matrix(bool (& e)[c_vvv][c_vv]) {
        for( int i = 0; i < c_v; ++i ) {
                for( int i2 = 0; i2 < c_v; ++i2 ) {
                        int const wi = (i * c_vv + i2 * c_v);
                        bool ek[c_v][c_vv] = {false};
                        get_Ek_matrix(ek, i * c_v);
                        for( int si = 0; si < c_v; ++si )
                                for( int sj = 0; sj < c_vv; ++sj )
                                        e[wi + si][sj] = ek[si][sj];
                }
        }
}

constexpr inline void get_Fk_matrix(bool (& fk)[c_v][c_vv], int const k) {
        for( int i = 0; i < c_v; ++i )
                for( int j = 0; j < c_vv; ++j )
                        fk[i][j] = false;
        for( int i = 0; i < c_v; ++i )
                fk[i][(i * c_v + k) % c_vv] = true;
}

constexpr inline void get_F_matrix(bool (& f)[c_vvv][c_vv]) {
        for( int i = 0; i < c_v; ++i ) {
                for( int i2 = 0; i2 < c_v; ++i2 ) {
                        int const wi = i * c_vv + i2 * c_v;
                        bool fk[c_v][c_vv] = {false};
                        get_Fk_matrix(fk, i2);

                        for( int si = 0; si < c_v; ++si )
                                for( int sj = 0; sj < c_vv; ++sj )
                                        e[wi + si][sj] = fk[si][sj];
                }
        }
}

constexpr inline void H2_matrix_write_at_row(bool (& h2)[3 * c_vv][c_vvv],
                                             bool const (& mt)[c_vvv][c_vv], int const row_idx) {
        for( int i = 0; i < c_vv; ++i ) {
                for( int j = 0; j < c_vvv; ++j ) {
                        h2[row_idx + i][j] = mt[j][i];
                }
        }
}

//ready to use BTW... but... they all are big :(
constexpr inline void get_H2_matrix(bool (& h2)[3 * c_vv][c_vvv]) {
        bool mt[c_vvv][c_vv] = {false};
        get_D_matrix(mt);
        H2_matrix_write_at_row(h2, mt, 0);
        get_E_matrix(mt);
        H2_matrix_write_at_row(h2, mt, c_v);
        get_F_matrix(mt);
        H2_matrix_write_at_row(h2, mt, cv + cv);
}


}
}

