#ifndef __itkMinimalPathImageFilter_txx
#define __itkMinimalPathImageFilter_txx

#include "itkMinimalPathImageFilter.h"
#include "itkConnectedComponentAlgorithm.h"
#include "itkConnectedComponentAlgorithmExtras.h"
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
  m_UseDistWeights = true;
  m_UnitCost = 0.00001;  // should be machine epsilon
  m_StartLabel = 1;
  m_EndLabel = 2;
  m_MarkLabel = 0;
  m_Direction = -1;
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

  //reg.SetSize(this->GetInput()->GetLargestPossibleRegion().GetSize());
  //reg.SetIndex(this->GetInput()->GetLargestPossibleRegion().GetIndex());
  //CostImage->SetRegions(reg);
  CostImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
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
    itkExceptionMacro(<< "Insufficient markers set - need at least start and end values");
    }
  if (m_MarkLabel == 0)
    {
    // not set by user
    m_MarkLabel = m_LabelChain[0];
    }

  // reset the cost vector
  m_CostVector.assign(0, 0.0);
  for (unsigned i=0;i<m_LabelChain.size()-1;i++)
    {
    if ((m_LabelChain[i] != 0) && (m_LabelChain[i+1] != 0))
      {
      ComputeLink(m_LabelChain[i], m_LabelChain[i+1], 
		  m_MarkLabel, CostImage, progress);
      }
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
  // loop over the label image
  bool FoundStart = false;
  bool FoundEnd   = false;
  IndexType StartIndex, EndIndex;

  while (!labIt.IsAtEnd())
    {
    // Find all pixels with value of the start label.
    LabelImagePixelType Val = labIt.Get();
    if (Val == StartLabel)
      {
      // point the other iterators to the correct spot
      IndexType Ind = labIt.GetIndex();
      if (!FoundStart)
	{
	FoundStart = true;
	// record the start index to compute the sign of desired direction
	StartIndex = Ind;
	}
      // mark the output label image
//       outIt.SetIndex(Ind);
//       outIt.Set(MarkLabel);
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
      if (!FoundEnd)
	{
	FoundEnd = true;
	EndIndex = Ind;
	}
//       outIt.SetIndex(Ind);
//       outIt.Set(MarkLabel);
      }
    ++labIt;
    }

  if (!FoundStart)
    {
    // exception - failed to find start label
    itkExceptionMacro(<< "Failed to find start label " << (int)StartLabel);
    }

  if (!FoundEnd)
    {
    itkExceptionMacro( << "Failed to find end label " << (int)EndLabel);
    }

  CNInputIterator inNIt(kernelRadius,
                        this->GetInput(),
                        this->GetOutput()->GetRequestedRegion() );


  // A label neighborhood iterator (only for initialization)
  typedef ConstShapedNeighborhoodIterator<LabelImageType> CNLabelIterator;
  CNLabelIterator labNIt(kernelRadius,
			 this->GetMarkerImage(),
			 this->GetOutput()->GetRequestedRegion() );

  typedef ShapedNeighborhoodIterator<LabelImageType> NLabelIterator;
  NLabelIterator outNIt(kernelRadius,
			this->GetOutput(),
			this->GetOutput()->GetRequestedRegion() );


  typedef ConstShapedNeighborhoodIterator<CostImageType> CNCostIterator;
  CNCostIterator costNIt(kernelRadius,
			 CostImage,
			 this->GetOutput()->GetRequestedRegion() );

  // figure out the direction information if necessary
  bool positive = true;

  if (m_Direction >= 0)
    {
    if (m_Direction >= TInputImage::ImageDimension)
      {
      itkExceptionMacro(<< "Illegal direction");
      }
    // The possible directions have been constrained
    positive = (EndIndex[m_Direction] - StartIndex[m_Direction]) > 0;
    setConnectivityDirection( &inNIt, m_Direction, positive);
    setConnectivityDirection( &labNIt, m_Direction, positive);
    setConnectivityDirection( &costNIt, m_Direction, positive);
    setConnectivityDirection( &outNIt, m_Direction, positive);
    }
  else
    {
    setConnectivity( &inNIt, m_FullyConnected );
    setConnectivity( &labNIt, m_FullyConnected );
    setConnectivity( &costNIt, m_FullyConnected );
    setConnectivity( &outNIt, m_FullyConnected );
    }

  ConstantBoundaryCondition<InputImageType> iBC;
  iBC.SetConstant(NumericTraits<InputImagePixelType>::max());
  inNIt.OverrideBoundaryCondition(&iBC);

  ConstantBoundaryCondition<LabelImageType> lBC;
  // This might seem like a strange boundary condition, but we'll
  // never be checking pixels with this label because they are already
  // minimum cost
  lBC.SetConstant(StartLabel);
  labNIt.OverrideBoundaryCondition(&lBC);
  outNIt.OverrideBoundaryCondition(&lBC);


  ConstantBoundaryCondition<CostImageType> cBC;
  cBC.SetConstant(NumericTraits<CostPixType>::max());
  costNIt.OverrideBoundaryCondition(&cBC);

  // set up weights
  WeightArrayType Weights;
  SetupWeights(Weights, inNIt);

  inNIt.GoToBegin();
  labNIt.GoToBegin();
  costNIt.GoToBegin();

  PixPriorityType TopPix;
  for (;;)
    {
    progress.CompletedPixel();
    // pop the head of the priority queue
    if (PriorityQueue.empty())
      {
      itkExceptionMacro(<< "Priority queue empty " << (int)StartLabel << " " << (int)EndLabel);
      }
    TopPix = PriorityQueue.top();
    PriorityQueue.pop();
    // set the iterators
    costIt.SetIndex(TopPix.location);
    labIt.SetIndex(TopPix.location);
    LabelImagePixelType ThisLab = labIt.Get();
    // check the label value - if it is EndValue, we are done.
    if (ThisLab == EndLabel) 
      {
      // record the cost 
      m_CostVector.push_back(TopPix.priority);
      break;
      }
    CostPixType ThisCost = costIt.Get();
    //std::cout << TopPix.location << " " << (int)ThisLab << " " << TopPix.priority << std::endl;
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
	   sIt != inNIt.End(); ++sIt, ++lIt, ++WIdx, ++cIt)
	{
	LabelImagePixelType NVal = lIt.Get();
	if (NVal != StartLabel)
	  {
	  InputImageOffsetType Off = sIt.GetNeighborhoodOffset();
	  PixPriorityType NewPix;
	  // the cost to reach the neighbor
	  NewPix.priority = (m_UnitCost+(CostPixType)sIt.Get()) * Weights[WIdx] + TopPix.priority;
	  if (NewPix.priority < cIt.Get())
	    {
	    // found a better path to the next pixel
	    NewPix.location = CentIndex + Off;
	    // insert into queue
	    PriorityQueue.push(NewPix);
	    }

	  }
	} 
      }
    }


//   typedef typename itk::ImageFileWriter<CostImageType> WType;
//   typename WType::Pointer wrt = WType::New();
//   wrt->SetInput(CostImage);
//   wrt->SetFileName("cost.nii");
//   wrt->Update();
  // backtrack
  if (m_Direction != -1)
    {
    // The possible directions have been constrained
    setConnectivityDirection( &inNIt, m_Direction, !positive);
    setConnectivityDirection( &labNIt, m_Direction, !positive);
    setConnectivityDirection( &costNIt, m_Direction, !positive);
    setConnectivityDirection( &outNIt, m_Direction, !positive);

    }

  IndexType CentIndex = TopPix.location;
  costNIt += CentIndex - costNIt.GetIndex();

#if 0
  // a method that is meant to deal with situations with lots of
  // identical/zero cost paths. Doesn't really work.
  labIt.SetIndex(CentIndex);
  outNIt += CentIndex - outNIt.GetIndex();
  outNIt.SetCenterPixel(MarkLabel);
  for (;;)
    {
      typename NLabelIterator::ConstIterator lIt;
      typename CNCostIterator::ConstIterator cIt;
      // look for the smallest neighbor in the cost image
      InputImageOffsetType MinOff;
      bool valid_neigh = false;
      // slightly more complex version which checks to see whether we
      // have already marked a voxel as being on the path to avoid
      // loops in zero cost zones
      CostPixType MN = NumericTraits<CostPixType>::max();
      //std::cout << CentIndex ;
      for (cIt = costNIt.Begin(); cIt != costNIt.End(); ++cIt)
	{
	// no good when getting to the last pixel
	  CostPixType NVal = cIt.Get();
	  //std::cout << NVal << " " << cIt.GetNeighborhoodOffset() << " " ;
	  if (NVal < MN)
	    {
	    valid_neigh = true;
	    MN = NVal;
	    MinOff = cIt.GetNeighborhoodOffset();
	    }
	}
      //std::cout << std::endl;
      if (!valid_neigh)
	{
	itkWarningMacro(<< "No valid neighbors - probably a zero cost region. Try increasing UnitCost");
	break;
	}
      CentIndex += MinOff;
      costNIt += MinOff;
      outNIt += MinOff;
      costNIt.SetCenterPixel(NumericTraits<CostPixType>::max());
      outNIt.SetCenterPixel(MarkLabel);
      labIt.SetIndex(CentIndex);
      progress.CompletedPixel();
      if (labIt.Get() == StartLabel)
	{
	// reached the start label
	break;
	}
    }
#else
  outIt.SetIndex(CentIndex);
  outIt.Set(MarkLabel);
  labIt.SetIndex(CentIndex);
  CostPixType MN = NumericTraits<CostPixType>::max();
  for (;;)
    {
      typename CNCostIterator::ConstIterator cIt;
      // look for the smallest neighbor in the cost image
      InputImageOffsetType MinOff;

      bool valid_neigh = false;
      for (cIt = costNIt.Begin(); cIt != costNIt.End(); ++cIt)
	{
	CostPixType NVal = cIt.Get();
	if (NVal < MN)
	  {
	  //std::cout << NVal << " " << MN << std::endl;
	  valid_neigh = true;
	  MN = NVal;
	  MinOff = cIt.GetNeighborhoodOffset();
	  }
	}
      if (!valid_neigh)
	{
	itkWarningMacro(<< "No valid neighbors - probably a zero cost region. Try increasing UnitCost");
	std::cout << MN << std::endl;
	break;
	}
      CentIndex += MinOff;
      costNIt += MinOff;
      outIt.SetIndex(CentIndex);
      outIt.Set(MarkLabel);
      labIt.SetIndex(CentIndex);
      progress.CompletedPixel();
      if (labIt.Get() == StartLabel)
	{
	// reached the start label
	break;
	}
    }
#endif
}

template <class TInputImage, class TLabelImage>
void
MinimalPathImageFilter<TInputImage, TLabelImage>
::SetupWeights(WeightArrayType &weights, const CNInputIterator &CNIt)
{
  typename CNInputIterator::ConstIterator sIt;
  InputImageSpacingType spacing = this->GetInput()->GetSpacing();
  WeightArrayType Tweights;
  for (sIt = CNIt.Begin(); !sIt.IsAtEnd(); ++sIt)
    {
    InputImageOffsetType Off = sIt.GetNeighborhoodOffset();
    CostPixType dist = 0;
    if (m_UseDistWeights)
      {
      for (int i = 0; i < ImageDimension; i++)
	{
	dist += (Off[i] * spacing[i]) * (Off[i] * spacing[i]);
	}
      dist = sqrt(dist);
      }
    else
      {
      dist = 1.0;
      }
    Tweights.push_back(dist);
//    std::cout << dist << " " << std::endl;
    }
  weights = Tweights;
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
