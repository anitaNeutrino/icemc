#include "icemc_random.h" 
#include "TUUID.h" 
#include "TRandom3.h" 

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
#include "TRandomGen.h" 
#endif

ClassImp(TRandomXoshiro256Plus); 



/**
 * A 64-bit state RNG 
 * from http://xoshiro.di.unimi.it/splitmix64.c 
 * Used to fully seed the total state
 */

static ULong_t next_ran_splitmix64(ULong_t x) 
{
    ULong_t z = (x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

/* see http://xoshiro.di.unimi.it/xoshiro256plus.c */ 
static inline ULong_t rotl(const ULong_t x, int k) 
{
  return (x << k) | (x >> (64-k)); 

}

/* Note that unlike in the reference implementation, we increment the state then return the new value
 * This is to avoid seeding of the state with e.g. run,event causing problems! 
 **/ 
ULong_t TRandomXoshiro256Plus::RawRndm() 
{
  const ULong_t t = fState[1] << 17; 
  fState[2] ^= fState[0];
  fState[3] ^= fState[1];
  fState[1] ^= fState[2];
  fState[0] ^= fState[3];
  fState[2] ^= t; 
  fState[3] = rotl(fState[3],45);

  return fState[0] + fState[3]; 
}

Double_t TRandomXoshiro256Plus::Rndm() 
{
  const Double_t inv_max = 1./(1ull << 53); 
  return (RawRndm() >> 11) * inv_max; 
}

void TRandomXoshiro256Plus::SetSeed(ULong_t seed) 
{
  if (seed == 0) 
  {
    TUUID uuid; 
    uuid.GetUUID( (UChar_t*) fState); 
    TUUID uuid2; 
    uuid2.GetUUID( ((UChar_t*) fState) + 16); 
  }
  else
  {
    fState[0] = next_ran_splitmix64(seed);
    fState[1] = next_ran_splitmix64(fState[0]); 
    fState[2] = next_ran_splitmix64(fState[1]); 
    fState[3] = next_ran_splitmix64(fState[2]); 
  }
}



void TRandomXoshiro256Plus::RndmArray(Int_t n, Float_t * arr) 
{
  const Float_t inv_max = 1./(1ul << 24); 
  for (Int_t i = 0; i < n/2; i++)
  {
    arr[2*i] = ((RawRndm() >> 16) & 0xffffff) *inv_max;
    arr[2*i+1] = (RawRndm() >>40)  *inv_max;
  }

  //we have one left over... 
  if (n % 2) arr[n-1] = (RawRndm() >>40)  *inv_max;

}

void TRandomXoshiro256Plus::RndmArray(Int_t n, Double_t * arr) 
{
  const Double_t inv_max = 1./(1ull << 53); 
  for (Int_t i = 0; i < n; i++) arr[i] = (RawRndm()>>11) *inv_max;
}




static ULong_t the_seed = 12345; 
static WhichIceMcRNGType the_rng_type; 
static TRandom* rngs[HOW_MANY_RNGS_DO_WE_HAVE]; 


class random_initializer
{
  public:
  random_initializer() 
  {
    setRNGType(RNG_TYPE_XOSHIRO256PLUS); 
  }
};

static random_initializer initializer; 


WhichIceMcRNGType getRNGType() { return the_rng_type; } 

void setRNGType(WhichIceMcRNGType type) 
{
  the_rng_type = type; 

  for (unsigned i = 0; i < HOW_MANY_RNGS_DO_WE_HAVE; i++) 
  {
    if (rngs[i]) delete rngs[i]; 

    rngs[i] = type == RNG_TYPE_XOSHIRO256PLUS ? (TRandom*) new TRandomXoshiro256Plus(the_seed) : 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,8)
              type == RNG_TYPE_MT64     ? (TRandom*) new TRandomMT64  : 
              type == RNG_TYPE_MIXMAX256     ? (TRandom*) new TRandomMixMax256  : 
#endif
              (TRandom*) new TRandom3; 
  }
}


void setSeed(ULong_t s) 
{
  the_seed = s; 
  for (unsigned i = 0; i <HOW_MANY_RNGS_DO_WE_HAVE; i++) rngs[i]->SetSeed(s + i * 1235); //each gets a different seed derived from the primary.
}


TRandom * getRNG(WhichIceMcRNG w) 
{
  return rngs[w]; 
}


