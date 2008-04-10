#ifndef __itkConnectedComponentAlgorithm_h
#define __itkConnectedComponentAlgorithm_h

#include "itkImage.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkShapedNeighborhoodIterator.h"

namespace itk
{

template< class TIterator >
TIterator*
setConnectivityDirection( TIterator* it, int direction, bool positive=true)
{
  // set a direction orthogonal to an axis - don't worry about fully
  // connected option
  typename TIterator::OffsetType offset;
  int val;
  it->ClearActiveList();
  if (positive)
    val = 1;
  else
    val = -1;

  unsigned int centerIndex = it->GetCenterNeighborhoodIndex();
  for( unsigned int d=0; d < centerIndex*2 + 1; d++ )
    {
    offset = it->GetOffset( d );
    if (offset[direction] == val) 
      {
      it->ActivateOffset( offset );
      }
    }
  return it;

}
}


#endif
