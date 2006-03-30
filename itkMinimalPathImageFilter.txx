#ifndef __itkMinimalPathImageFilter_txx
#define __itkMinimalPathImageFilter_txx

#include "itkMinimalPathImageFilter.h"
#include "itkConnectedComponentAlgorithm.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkProgressReporter.h"
#include "itkImageFileWriter.h"

namespace itk {

template <class TInputImage, class TLabelImage>
MinimalPathImageFilter<TInputImage, TLabelImage>
::MinimalPathImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  m_FullyConnected = true;
  m_UnitCost = 1.0;
  m_StartLabel = 1;
  m_EndLabel = 2;
  m_MarkLabel = 0;
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
  // set up the cost image
  typename CostImageType::Pointer CostImage = CostImageType::New();
  typename CostImageType::RegionType reg;

  reg.SetSize(this->GetInput()->GetLargestPossibleRegion().GetSize());
  reg.SetIndex(this->GetInput()->GetLargestPossibleRegion().GetIndex());
  CostImage->SetRegions(reg);
  CostImage->Allocate();
  this->AllocateOutputs();

  this->GetOutput()->FillBuffer(0);

  ProgressReporter progress(this, 0, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels()*2);

  // set up the labels defining the path
  if (m_LabelChain.size() == 0)
    {
    // use StartLabel and EndLabel
    m_LabelChain.push_back(m_StartLabel);
    m_LabelChain.push_back(m_EndLabel);
    }

  if (m_LabelChain.size() < 2)
    {
    // exception - not enough labels set
    }
  if (m_MarkLabel == 0)
    {
    // not set by user
    m_MarkLabel = m_LabelChain[0];
    }

  for (int i=0;i<m_LabelChain.size()-1;i++)
    {
    ComputeLink(m_LabelChain[i], m_LabelChain[i+1], 
		m_MarkLabel, CostImage, progress);
    }
}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::ComputeLink(const LabelImagePixelType StartLabel, 
	      const LabelImagePixelType EndLabel, 
	      const LabelImagePixelType MarkLabel,
	      typename CostImageType::Pointer CostImage,
	      ProgressReporter &progress)
{
  // compute a link in our path
  // initialize the cost image (allocated in calling procedure
  CostImage->FillBuffer(NumericTraits<CostPixType>::max());

  // prepare iterators
  typedef ImageRegionConstIterator<LabelImageType> LabIterator;
  typedef ImageRegionIterator< CostImageType> CostIterator;
  typedef ImageRegionIterator< LabelImageType> OutIterator;
  typedef ImageRegionConstIterator< InputImageType> InputIterator;
    // neighborhood iterator
  InputImageSizeType kernelRadius;
  kernelRadius.Fill(1);

  CNInputIterator inNIt(kernelRadius,
                        this->GetInput(),
                        this->GetOutput()->GetRequestedRegion() );

  setConnectivity( &inNIt, m_FullyConnected );
  ConstantBoundaryCondition<InputImageType> iBC;
  iBC.SetConstant(NumericTraits<InputImagePixelType>::max());
  inNIt.OverrideBoundaryCondition(&iBC);
  // set up weights
  WeightArrayType Weights;
  SetupWeights(Weights, inNIt);
  // A label neighborhood iterator (only for initialization)
  typedef ConstShapedNeighborhoodIterator<LabelImageType> CNLabelIterator;
  CNLabelIterator labNIt(kernelRadius,
			 this->GetMarkerImage(),
			 this->GetOutput()->GetRequestedRegion() );

  setConnectivity( &labNIt, m_FullyConnected );
  ConstantBoundaryCondition<LabelImageType> lBC;
  // This might seem like a strange boundary condition, but we'll
  // never be checking pixels with this label because they are already
  // minimum cost
  lBC.SetConstant(StartLabel);
  labNIt.OverrideBoundaryCondition(&lBC);

  typedef ConstShapedNeighborhoodIterator<CostImageType> CNCostIterator;
  CNCostIterator costNIt(kernelRadius,
			 CostImage,
			 this->GetOutput()->GetRequestedRegion() );
  setConnectivity( &costNIt, m_FullyConnected );
  ConstantBoundaryCondition<CostImageType> cBC;
  cBC.SetConstant(NumericTraits<CostPixType>::max());
  costNIt.OverrideBoundaryCondition(&cBC);

  // create the queue
  PriorityQueueType PriorityQueue;
  // initialise the queue
  LabIterator labIt( this->GetMarkerImage(),
		     this->GetOutput()->GetRequestedRegion() );
  OutIterator outIt( this->GetOutput(),
		     this->GetOutput()->GetRequestedRegion() );
  CostIterator costIt( CostImage,
		       this->GetOutput()->GetRequestedRegion() );
  InputIterator inIt(this->GetInput(),
		     this->GetOutput()->GetRequestedRegion() );
  
  labIt.GoToBegin();
  costIt.GoToBegin();
  outIt.GoToBegin();
  inIt.GoToBegin();
  inNIt.GoToBegin();
  labNIt.GoToBegin();
  costNIt.GoToBegin();
  // loop over the label image
  bool FoundStart = false;
  bool FoundEnd   = false;
  while (!labIt.IsAtEnd())
    {
    // Find all pixels with value of the start label.
    LabelImagePixelType Val = labIt.Get();
    if (Val == StartLabel)
      {
      // point the other iterators to the correct spot
      FoundStart = true;
      IndexType Ind = labIt.GetIndex();
      inIt.SetIndex(Ind);
      // mark the output label image
      outIt.SetIndex(Ind);
      outIt.Set(MarkLabel);
      // Insert in the queue (equal priority, 0)
      // Don't worry about finding pixels on the edge of the marker at
      // this stage because those issues will be dealt with when
      // processing the queue.
      PixPriorityType ThisPix;
      ThisPix.location=Ind;
      ThisPix.priority=0;
      PriorityQueue.push(ThisPix);
      } 
    else if (Val == EndLabel)
      {
      IndexType Ind = labIt.GetIndex();
      outIt.SetIndex(Ind);
      outIt.Set(MarkLabel);
      FoundEnd = true;
      }
    ++labIt;
    }

  if (!FoundStart)
    {
    // exception - failed to find start label
    }

  if (!FoundEnd)
    {
    // exception - failed to find end label
    }

  PixPriorityType TopPix;
  for (;;)
    {
    progress.CompletedPixel();
    // pop the head of the priority queue
    TopPix = PriorityQueue.top();
    PriorityQueue.pop();
    // set the iterators
    costIt.SetIndex(TopPix.location);
    labIt.SetIndex(TopPix.location);
    LabelImagePixelType ThisLab = labIt.Get();
    // check the label value - if it is EndValue, we are done.
    if (ThisLab == EndLabel) break;

    CostPixType ThisCost = costIt.Get();
    // std::cout << TopPix.location << " " << (int)ThisLab << " " << TopPix.priority << std::endl;
    if (TopPix.priority < ThisCost)
      {
      labNIt += TopPix.location - labNIt.GetIndex();
      inNIt += TopPix.location - inNIt.GetIndex();
      costNIt += TopPix.location - costNIt.GetIndex();
      // if cost on queue is less than current cost of reaching that
      // location reset the cost of reacing that location.
      costIt.Set(TopPix.priority);
      // Evaluate the costs of reaching neighbors and insert in priority
      // queue
      typename CNInputIterator::ConstIterator sIt;
      typename CNLabelIterator::ConstIterator lIt;
      typename CNCostIterator::ConstIterator cIt;
      int WIdx;
      IndexType CentIndex = inNIt.GetIndex();
      for (WIdx=0,sIt = inNIt.Begin(), lIt = labNIt.Begin(), cIt = costNIt.Begin();
	   !sIt.IsAtEnd(); ++sIt, ++lIt, ++WIdx, ++cIt)
	{
	LabelImagePixelType NVal = lIt.Get();
	if (NVal != StartLabel)
	  {
	  InputImageOffsetType Off = sIt.GetNeighborhoodOffset();
	  PixPriorityType NewPix;
	  // the cost to reach the neighbor
	  NewPix.priority = (m_UnitCost+sIt.Get()) * Weights[WIdx] + TopPix.priority;
	  if (NewPix.priority < cIt.Get())
	    {
	    // found a better path to the next pixel
	    NewPix.location = CentIndex + Off;
	    // insert into queue
	    PriorityQueue.push(NewPix);
	    //std::cout << "pushing " << NewPix.location << " " << NewPix.priority << std::endl;
	    }

	  }
	} 
      }
    }

  // backtrack
  IndexType CentIndex = TopPix.location;
  costNIt += CentIndex - costNIt.GetIndex();
  outIt.SetIndex(CentIndex);
  outIt.Set(MarkLabel);
  for (;;)
    {
      typename CNCostIterator::ConstIterator cIt;
      // look for the smallest neighbor in the cost image
      CostPixType MN = NumericTraits<CostPixType>::max();
      InputImageOffsetType MinOff;
      for (cIt = costNIt.Begin();
	   !cIt.IsAtEnd(); ++cIt)
	{
	CostPixType NVal = cIt.Get();
	//std::cout << NVal << " " << cIt.GetNeighborhoodOffset() << " " ;
	if (NVal < MN)
	  {
	  MN = NVal;
	  MinOff = cIt.GetNeighborhoodOffset();
	  }
	}
      //std::cout << std::endl;
      //std::cout << "Backtrack " << CentIndex << " " << MN << std::endl;
      CentIndex += MinOff;
      costNIt += MinOff;
      outIt.SetIndex(CentIndex);
      outIt.Set(MarkLabel);
      progress.CompletedPixel();
      if (MN == 0)
	{
	// reached the start label
	break;
	}
    }

}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::SetupWeights(WeightArrayType &weights, const CNInputIterator &CNIt)
{
  typename CNInputIterator::ConstIterator sIt;
  InputImageSpacingType spacing = this->GetInput()->GetSpacing();
  for (sIt = CNIt.Begin(); !sIt.IsAtEnd(); ++sIt)
    {
    InputImageOffsetType Off = sIt.GetNeighborhoodOffset();
    CostPixType dist = 0;
    for (int i = 0; i < ImageDimension; i++)
      {
      dist += (Off[i] * spacing[i]) * (Off[i] * spacing[i]);
      }
    dist = sqrt(dist);
    weights.push_back(dist);
//    std::cout << dist << " " << std::endl;
    }
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
