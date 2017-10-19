# ldpc
LDPC (Low Density Parity Check) codes as implemented by me. both soft &amp; hard decoding. C++14 and above; Header-only library.

# What is this for ?
They allow to recover partially damaged information. They are especially good if you roughly know probability that some bit was damaged during a transfer (for example by measuring signal level).
I found there is really no libraries free to use with when LDPC codeword is not too long. Everyone it seems tries to use as long as possible codes in order to raise their recover abilities, but for my task I needed small ones. Also I found it's very hard read corresponding articles in the topic, all the scientiest here speak their bad understood by regular developers language... So I had to learn it a little in order to write the implementation and understand others.

# Theory
LDPC codes work on per-bit basis. I.e. from viewpoint of code there are no bytes, just series of bits.
In order to make encoding/decoding you need two matrices: decoding matrix, usually it's called H and encoding G, or Generator Matrix. When you build a code you put some number of bits to be an information for recovery. For matrix H the number of columns is equal to total codeword length (number of bits you are going to transfer via noisy channel), the number of rows is equal to number of check bits (i.e. these exceess bits for recovery), the length of useful message (or what you are going to transfer) then is equal <number of columns> minus <number of rows>
In the code H matrix is stored in sparsed form, which means that a number in a row shows a column where bit is equal to 1, -1 indicates end of the signature.
  
  There are three kind of decoding supported by LDPC:
  BEC, Hard decision and Soft. BEC means binary erasure channel, when bit cannot change its value from 1 to 0 or from 0 to 1, but you can know that this concrete bit was unrecognized on receiving. I.e. you dont know its initial value. This algorithm is not implemented in the library because I didn't need it.
  Hard decision is most simple, and Soft allows you to give probabilities that some bit is was recognized in a right way.
  
  'Soft' decoding IMHO works much better. You should provide 'float' probabilities that the bit is equal to 1. I.e. if it's 1 then the bit is recognized for sure correctly and an errror here is impossible. If it's 0 then the bit is 0 without any doubts. 0.5 means that you have no idea what original bit was. 0.7 means that likely bit was 1 but you are unsure.
  
# How to use ?
There is a little test written, so it should help you somehow.
You need only ldpc folder, the library is header only
I hope I'll write in details soon.

# Where to take matrices ?
You can generate any with help of utilities [from here](http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/pchk.html) or by using my ldpc_build_matrix_v02.hpp it's advantage is that it can eliminate loops of any length on the expense of matrix computation time duration. However my matrix builder cannot create inverse generator (G) matrices so you either should write your own implementation or again use make-gen tool from this http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/progs.html 
+ regular expressions to make the matrix embeddable in C++ programm.
Also you can borrow matrices from places like https://github.com/tracierenea/GNU-Radio-GSoC2013/tree/master/alist_files or other collections. FWIW I currently use my own generated with ldpc_build_matrix_v02.hpp
In ldpc folder there is a small matrix collection except ldpc_matrix.hpp, there is ldpc_matrix.300x152x3x6.hpp and ldpc_matrix.n116.k64.wc3.l6_0.hpp

# License
MIT - or free to use in any form including commercial products. The author would appreciate if you mention him or his library, but it's not mandatory and I will not be offended in any way

# TODO
- In theory the library can be optimized in times if use bit logic instead of booleans.
- There is no way to compute G-matrices... I just had not mind force to writ it.
- BEC (Binary erasure channel) algorithm
