clear all
close all

syms s1 s2 s3 s4;

S=[s1 -s2 -s3 -s4 conj(s1) -conj(s2) -conj(s3) -conj(s4);
   s2  s1  s4 -s3 conj(s2)  conj(s1)  conj(s4) -conj(s3);
   s3 -s4  s1  s2 conj(s3) -conj(s4)  conj(s1)  conj(s2)]

E=S*S'

