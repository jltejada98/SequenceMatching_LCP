/*
 * SuffixArray.hpp
 * 
 * Copyright 2015 Robert Bakaric <rbakaric@irb.hr>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include <vector>
#include <string>
#include "Sais.h"

using namespace std;



/* SA container */

template <typename Tint>
class SuffixArray{

  protected:
/* Data Container */
  vector<Tint> SA;
  

/* Functions */ 
  void ComputeSuffixArray(const string& text);

  public:

/* Constructor */
  SuffixArray(const string & text);
/* Destructor */
  ~SuffixArray();
/* Explicite Constructor */
  void make(const string & text);
/* Explicite destructor */
  void destroy();
/* Getters */
  vector<Tint>& GetSuffixArray();
};



/* Constructors */

template <typename Tint>
SuffixArray<Tint>::SuffixArray(const string & text){
  ComputeSuffixArray(text);
}

template <typename Tint>
void SuffixArray<Tint>::make(const string & text){
  ComputeSuffixArray(text);
}



/* Destructors */

template <typename Tint>
SuffixArray<Tint>::~SuffixArray(){
  SA.clear();
}

template <typename Tint>
void SuffixArray<Tint>::destroy(){
  SA.clear();
}



/* SA computation */

template <typename Tint>
void SuffixArray<Tint>::ComputeSuffixArray(const string& text){
//   std::cerr << "Text size:" << text.size() << std::endl;
   SA.resize(static_cast<Tint>(text.size()));
   vector<unsigned char> data(text.begin(), text.end());
   if(sais(&data[0], &SA[0], static_cast<Tint>(text.size())) != 0)
      throw runtime_error ("Cannot allocate enough memory \n" );
}


/* Getter */
template <typename Tint>
vector<Tint> & SuffixArray<Tint>::GetSuffixArray(){
   
   return SA;
}

