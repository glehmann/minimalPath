#ifndef __itkMinimalPathImageFilter_txx
#define __itkMinimalPathImageFilter_txx

#include "itkMinimalPathImageFilter.h"

namespace itk {

template <class TInputImage, class TLabelImage>
MinimalPathImageFilter<TInputImage, TLabelImage>
::MinimalPathImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  m_FullyConnected = true;
}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  LabelImagePointer  markerPtr = this->GetMarkerImage();

  InputImagePointer  inputPtr =
    const_cast< InputImageType * >( this->GetInput() );

  if ( !markerPtr || !inputPtr )
    {
    return;
    }
  // We need to
  // configure the inputs such that all the data is available.
  //
  markerPtr->SetRequestedRegion(markerPtr->GetLargestPossibleRegion());
  inputPtr->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()
    ->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::GenerateData()
{
  this->AllocateOutputs();
}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
}

// end namespace itk
}

#endif
