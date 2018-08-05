#ifndef ICEMC_LOCAL_COORDINATES_H
#define ICEMC_LOCAL_COORDINATES_H


#include "Geoid.h"
#include "TVector3.h"

namespace icemc {

  /**
   * A class for constructing local coordinate systems centered on a Geoid::Position
   * The z-axis will always point directly upwards from the origin position.
   */

  class LocalCoordinateSystem {
  public:
    
    LocalCoordinateSystem(const Geoid::Position& origin)
      : fOrigin(origin){
      init();
    }
    LocalCoordinateSystem(const Geoid::Position& origin,  const Geoid::Position& suggestedX)
      : fOrigin(origin){
      init(&suggestedX);
    }

    LocalCoordinateSystem(const Geoid::Position& origin,  const Geoid::Position& suggestedX,  const TVector3& newZ)
      : fOrigin(origin), fLocalZ(newZ.Unit()){
      init(&suggestedX);
    }
    
    inline Geoid::Position localPositionToGlobal(const TVector3& localPos) const {
      return toGlobal(localPos, true);
    }

    inline TVector3 localTranslationToGlobal(const TVector3& localTranslation) const {
      return toGlobal(localTranslation,  false);
    }
    
    inline TVector3 globalPositionToLocal(const Geoid::Position& globalPos) const {
      return toLocal(globalPos, true);
    }

    inline TVector3 globalTranslationToLocal(const TVector3& globalTranslation) const {
      return toLocal(globalTranslation,  false);
    }    
    

  private:

    const Geoid::Position fOrigin;
    TVector3 fLocalX; // unit vector describing local x-axis in global coordinates
    TVector3 fLocalY; // unit vector describing local y-axis in global coordinates
    TVector3 fLocalZ; // unit vector describing local z-axis in global coordinates
    


    /**
     * Set up the local coordinate system basis unit vectors
     */
    inline void init(const Geoid::Position* optionalX = NULL){

      if(fLocalZ.Mag()==0){
	Geoid::Position up = fOrigin;
	up.SetAltitude(fOrigin.Altitude()+1);
	fLocalZ = (up - fOrigin).Unit(); // directly above local origin
      }

      if(optionalX){
	fLocalY = fLocalZ.Cross(*optionalX - fOrigin).Unit(); // y is perpendicual to local up AND the direction of optionalX
      }
      else{
	fLocalY = fLocalZ.Orthogonal().Unit(); // any old perpendicular vector
      }
      fLocalX = fLocalY.Cross(fLocalZ).Unit(); // ... and therefore x is local horizontal in direction of optionalX, if provided
    }


    /** 
     * Workhorse transformation from global to local
     * 
     * @param v the vector-ish thing to transform
     * @param translate do true for positions, do false for translations
     * 
     * @return v in the local coordinate system
     */
    inline TVector3 toLocal(const TVector3& v, bool translate) const {
      // translate the earth centered vector, v, to our new local coordinate system
      int translationFactor = translate ? 1 : 0;
      const TVector3 v2 = v - translationFactor*fOrigin;
      double x = v2.Dot(fLocalX);
      double y = v2.Dot(fLocalY);
      double z = v2.Dot(fLocalZ);
      return TVector3(x, y, z);
    }

    /** 
     * Workhorse transformation from local to global
     * 
     * @param v the vector-ish thing to transform
     * @param translate do true for positions, do false for translations
     * 
     * @return v in the global coordinate system
     */
    inline TVector3 toGlobal(const TVector3& v, bool translate) const {
      TVector3 v2(0, 0, 0);
      int translationFactor = translate ? 1 : 0;
      v2 += translationFactor*fOrigin;
      v2 += v.X()*fLocalX;
      v2 += v.Y()*fLocalY;
      v2 += v.Z()*fLocalZ;
      return v2;
    }
    

  };
  
  

};



#endif // ICEMC_LOCAL_COORDINATES_H
