/** \file
 *
 *  $Date: $
 *  $Revision: $
 */

#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"

using namespace GeomDetEnumerators;

std::ostream& operator<<( std::ostream& s, Location l) {
  if (l == barrel) return s << "barrel";
  else return s << "endcap";
}

std::ostream& operator<<( std::ostream& s, SubDetector m){
  if ( m == PixelBarrel ) return s << "PixelBarrel";
  else if ( m == PixelEndcap ) return s << "PixelEndcap";
  else if ( m == TIB) return s << "TIB";
  else if (m == TOB) return s << "TOB";
  else if (m == TID) return s << "TID";
  else if (m == TEC) return s << "TEC";
  else if ( m == DT ) return s << "DT";
  else if ( m == CSC ) return s << "CSC";
  else if ( m == RPCBarrel ) return s << "RPCBarrel";
  else if ( m == RPCEndcap ) return s << "RPCEndcap";
  else return s << "?";
}


