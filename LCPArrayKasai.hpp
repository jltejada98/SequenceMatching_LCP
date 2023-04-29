/*
 * LCPArrayKasai.hpp
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
 * 
 */

#include <vector>

template <typename Tint>
class LCPArrayKasai {
  
   vector <Tint> LCP;
   vector <Tint> Rank;

  
   void ComputeLCPArray(const vector<Tint>& SA, string& text);
   public:
/* Constructor */
   LCPArrayKasai(const vector<Tint>& SA, string& text);
   
/* Explicit Constructor */
   void make(const vector<Tint>& SA, string& text);
   
/* Destructor */
   ~LCPArrayKasai();

/* Explicit  Destructor */
   void destroy();
   
   vector<Tint> GetLCPArray();
};


template<typename Tint>
LCPArrayKasai<Tint>::LCPArrayKasai(const vector<Tint>& SA, string& text){
   ComputeLCPArray(SA, text);
}

template<typename Tint>
void LCPArrayKasai<Tint>::make(const vector<Tint>& SA, string& text){
   ComputeLCPArray(SA, text);
}

template<typename Tint>
LCPArrayKasai<Tint>::~LCPArrayKasai(){
   LCP.clear();
   Rank.clear();
}

template<typename Tint>
void LCPArrayKasai<Tint>::destroy(){
   LCP.clear();
   Rank.clear();
}


/* Functions */
template<typename Tint>
void LCPArrayKasai<Tint>::ComputeLCPArray(const vector<Tint>& SA, string& text){
 
   /*Kasai et al. 2001.*/
   Rank.resize(SA.size());
   LCP.resize(SA.size());
   
   for (Tint i =0; i< SA.size(); i++)
      Rank[SA[i]] = i;
      
   Tint lcp=0;
  
   for(Tint i=0; i< SA.size(); i++){
      if( Rank[i] > 0){
         Tint j = SA[Rank[i]-1];
         while(text[i+lcp] == text[j+lcp])
            lcp++;
         LCP[Rank[i]] = lcp;
         if(lcp>0) lcp--;
      }
   }
}

template <typename Tint>
vector<Tint> LCPArrayKasai<Tint>::GetLCPArray(){
   return LCP;
} 
