#pragma once

namespace ldpc {

//in most compact form
template<int m, int n>
void print_matrix(bool const (& matrix)[m][n]) {
        printf("%02dx%02d binary matrix: \n", m, n);
        for( int i = 0; i < m; ++i ) {
                for( int j = 0; j < n; ++j )
                        putc(matrix[i][j] ? '1' : '0', stdout);
                putc('\n', stdout);
        }
}

template <int m, int n>
void print_matrix_c_sparsed(bool const (&mx)[m][n])
{
        int wr = 0;
        for( int i = 0; i < m; ++i ) {
                int cur_wr = 0;
                for( int j = 0; j < n; ++j )
                        cur_wr += mx[i][j];
                if (cur_wr > wr) wr = cur_wr;
        }

        printf("C-style sparsed matrix %dx%d max row weight including -1: %d\n", m, n, wr + 1);
        for( int i = 0; i < m; ++i ) {
                printf("{");
                for( int j = 0; j < n; ++j ) {
                        if (mx[i][j]) printf("%d, ", j);
                }
                printf("-1}%s\n", (i+1 != m)? "," : "");
        }
};

//matrix in form suitable for make-pchk utility
//you can take it here: http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/pchk.html
template <int m, int n>
void print_matrix_for_make_pchk(bool const (&mx)[m][n]) {
        printf("matrix for make-pchk %dx%d\n", m, n);
        printf("make-pchk matrix.h.mx %d %d ", m, n);
        for( int i = 0; i < m; ++i )
                for( int j = 0; j < n; ++j)
                        if (mx[i][j])
                                printf("%d:%d ", i, j);
}


}