#pragma once
#include <stdio.h>
#include <math.h>
#include <vector>

namespace bvv
{
template <class T, int BUFFER_DEPTH = 2>
/*
  Class TBuffer.
  - Create class:
    : TBuffer <double, 9> b(1.23); # Buffer depth is 9.
    : TBuffer <double> b(1.23);    # Default buffer depth is 2.
    : TBuffer <double> b = 1.23; // DISCOURAGED: The behaviour may differ with different compilers!
    : TBuffer <double> b; # <-- defaults to zero.
    In all cases b obtains a status of "NOT modified", see next entry.
  - Check if modified:
    : if (b.modified) blah-blah
    : if (*b) blah-blah
    This functionality was the primary motivation for creation of the class.
    It enables you to keep an old value and the most recent one obtained from the last assignment
    in one place, under the same name, say "a", "b" and "c" instead of "a_new, a_old", etc. and do:
    : if (*a || *b || *c) if anything has changed in the controlled set of variables 
  - Keep history of "unique" assignments to the object: 
    std::vector<double> history = b.history();
    - "Unique" means: assignment sequence "x x y x" will be remembered as x y x.
 
*/
class TBuffer {
    T buffer[BUFFER_DEPTH];
    int history_len;
    int ind;
  public:
    // static const int BUF_DEPTH = 2;
    bool modified;

    int max_size() {return BUFFER_DEPTH;}

    TBuffer (T val=0)
    {
       ind = 0;
       buffer[0] = val;
       history_len = 1;
       modified = false;
    }

    T get_val()
    {
      return buffer[ind]; 
    }

    int get_hist_len()
    {
      return history_len; 
    }

    void set_val(T val)
    {
      T current_val = buffer[ind]; 
      if (current_val == val) {
          modified = false;
          return;
      }
      ind = (ind + 1) % BUFFER_DEPTH;
      // Buffer can only contain BUFFER_DEPTH elements:
      if (history_len < BUFFER_DEPTH) history_len++; 
      buffer[ind] = val;
      modified = true;
    }

    int prev_ind(){
      int ind_minus_one = ind - 1;
      if (ind_minus_one < 0) ind_minus_one = BUFFER_DEPTH - 1;
      return ind_minus_one;
    }

    std::vector<T> history(){
      std::vector<T> history(history_len);
      int buf_ind = ind;
      for (int i = 0; i < history_len; i++) {
        history[i] = buffer[buf_ind--];
        // Wrapping around:
        if (buf_ind < 0) buf_ind = BUFFER_DEPTH - 1;
      }
      return history;
    }

//    explicit operator bool() {
//      printf("changing!\n");
//      return modified;
//     }

    bool operator*() {
      return modified;
    }

    TBuffer<T, BUFFER_DEPTH>& operator=(const T& other) // copy assignment
    {
      set_val(other);
      return *this;
    }

   // To do things like double(x) or int y = x where x is of type TBuffer type.
   operator T() {return get_val();}
};
}
