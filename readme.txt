How to use dragon to get dynamic information about call graph and control
flow graph.

1. use opencc/openf90 to compile your program with feedback option
 opencc -fb_create myfeedback -fb_type=1 -fb_phase=0 test.c

2. run program using some  input 
 a.out
 
3. dump call graph/control flow graph using -dragon option
  opencc -fb_opt myfeedback -dragon -ipa -O2 test.c

4. invoke dragon to show graph

To get static information

Generating call graph and control flow graph
1. opencc -dragon -ipa -O2 test.c

Generating dependence information
2. opencc -dragon -O3 test.c


By Chunhua Liao ,4/1/2004
   liaoch@cs.uh.edu



