#ifndef __itkMinimalPathImageFilter_h
#define __itkMinimalPathImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkProgressReporter.h"

#include <queue>
#include <vector>

namespace itk {

/** \class MinimalPathImageFilter
 * \brief Find the lowest cost path through an image.
 *
 * Locations along the path are defined using a marker
 * image. Successive links in the path are computed between labels
 * defined in either StartLabel and EndLabel or in the vector
 * LabelChain. The paths are computed independently, and could
 * cross. The path is marked in the output image with the label
 * corresponding to the start of the path. Pixels in the input image
 * are used to compute costs, so they will need to be non
 * negative. The spacing information in the input image is used to
 * weight costs.
 *
 * \author Richard Beare, Department of Medicine, Monash University,
 * Melbourne, Australia
 *
 */
template<class TInputImage, class TLabelImage>
class ITK_EXPORT MinimalPathImageFilter : 
    public ImageToImageFilter<TInputImage, TLabelImage>
{
public:
  /** Standard class typedefs. */
  typedef MinimalPathImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TLabelImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef TLabelImage LabelImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::OffsetType     InputImageOffsetType;
  typedef typename InputImageType::SpacingType    InputImageSpacingType;
  typedef typename LabelImageType::Pointer        LabelImagePointer;
  typedef typename LabelImageType::ConstPointer   LabelImageConstPointer;
  typedef typename LabelImageType::RegionType     LabelImageRegionType;
  typedef typename LabelImageType::PixelType      LabelImagePixelType;
  
  typedef typename LabelImageType::IndexType      IndexType;
  /** Used in the internal buffer to accumulate costs */
  typedef double CostPixType;

  /** A type to support a set of labels, rather than just start and
   * end */
  typedef typename std::vector<LabelImagePixelType> LabelVectorType;
  typedef typename std::vector<CostPixType> CostVectorType;

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(MinimalPathImageFilter, 
               ImageToImageFilter);
  

   /** Set the marker image */
  void SetMarkerImage(TLabelImage *input)
     {
     // Process object is not const-correct so the const casting is required.
     this->SetNthInput( 1, const_cast<TLabelImage *>(input) );
     }

  /** Get the marker image */
  LabelImageType * GetMarkerImage()
    {
    return static_cast<LabelImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

   /** Set the input image */
  void SetInput1(TInputImage *input)
     {
     this->SetInput( input );
     }

   /** Set the marker image */
  void SetInput2(TLabelImage *input)
     {
     this->SetMarkerImage( input );
     }

  /**
   * Set/Get whether the path is defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOn. This only make sense if the direction isn't
   * constrained, and isn't particularly useful even then
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);

  /**
   * Default = true. Should costs be weighted according to the spacing information -
   * if yes then it will cost more to cross a diagonal, and the costs
   * in different directions will be different if the voxels aren't
   * cubic
  */
  itkSetMacro(UseDistWeights, bool);
  itkGetConstReferenceMacro(UseDistWeights, bool);
  itkBooleanMacro(UseDistWeights);

  itkSetMacro(StartLabel, LabelImagePixelType);
  itkGetConstReferenceMacro(StartLabel, LabelImagePixelType);

  itkSetMacro(EndLabel, LabelImagePixelType);
  itkGetConstReferenceMacro(EndLabel, LabelImagePixelType);

  /** If this isn't set the first element of LabelChain is used
   * instead */
  itkSetMacro(MarkLabel, LabelImagePixelType);
  itkGetConstReferenceMacro(MarkLabel, LabelImagePixelType);

  void SetLabelChain(const LabelVectorType &LV)
  {
    m_LabelChain = LV;
  }

  const LabelVectorType & GetLabelChain() const
  {
    return m_LabelChain;
  }

  const CostVectorType & GetCosts() const
  {
    return m_CostVector;
  }


  /** 
   * The UnitCost defines the cost of moving over pixels with
   * a brightness of zero. If there are lots of zero cost pixels then
   * the backtracking becomes difficult because it is easy to loop
   * back on yourself in the flat zone. It is therefore best to
   * provide a low cost of crossing a pixel. This can be done by
   * adding a constant to the image, but this method lets you do it
   * without changing the type of the input image.
   * The default value is  0.00001. */
  itkSetMacro(UnitCost, CostPixType);
  itkGetConstReferenceMacro(UnitCost, CostPixType);

  itkSetMacro(Direction, int);
  itkGetConstReferenceMacro(Direction, int);

protected:
  MinimalPathImageFilter();
  ~MinimalPathImageFilter() {};
  
  // the image that is going to be used to accumulate the cost
  // Could be a template parameter, but there doesn't appear to be
  // much point. Easiest to use a floating point type since there is
  // going to be weighting by chamfer distances. Otherwise it would be
  // necessary to mess around with scaling.
  typedef typename itk::Image<CostPixType, TInputImage::ImageDimension> CostImageType;

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** MinimalPathImageFilter needs  the entire input be
   * available. Thus, it needs to provide an implementation of
   * GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** This filter will enlarge the output requested region to produce
   * all of the output if the filter is configured to run to
   * convergence.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

  /** Single-threaded version of GenerateData. */
  void GenerateData();
  


private:
  MinimalPathImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_FullyConnected, m_UseDistWeights;
  LabelImagePixelType m_StartLabel, m_EndLabel, m_MarkLabel;
  LabelVectorType m_LabelChain;
  CostVectorType m_CostVector;
  CostPixType m_UnitCost;
  int m_Direction;

  typedef ConstShapedNeighborhoodIterator<InputImageType> CNInputIterator;
  typedef std::vector<CostPixType> WeightArrayType;

  void SetupWeights(WeightArrayType &weights, const CNInputIterator &CNIt);

  void ComputeLink(const LabelImagePixelType StartLabel, 
		   const LabelImagePixelType EndLabel, 
		   const LabelImagePixelType MarkLabel,
		   typename CostImageType::Pointer CostImage,
		   ProgressReporter &progress);

  // support for the priority queue
  typedef class PixPriorityType 
  {
  public:
    IndexType location;
    CostPixType priority;
  };

  class PixPriorityCompare
  {
  public:
    bool operator()(PixPriorityType const &l,
		    PixPriorityType const &r)
    {
      return l.priority > r.priority;
    }
  };

  typedef std::priority_queue<PixPriorityType, std::vector<PixPriorityType>, PixPriorityCompare > PriorityQueueType;
//  typedef std::priority_queue<PixPriorityType, std::deque<PixPriorityType>, PixPriorityCompare > PriorityQueueType;

} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMinimalPathImageFilter.txx"
#endif

#endif


