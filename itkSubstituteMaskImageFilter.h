#ifndef __itkSubstituteMaskImageFilter_h
#define __itkSubstituteMaskImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"


namespace itk
{
  
/** \class SubstituteMaskImageFilter
 * \brief Each location that is nonzero in a corresponding mask image
 * is replaced by a user defined value. Used for reseting regions
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TInput1, class TInput2=TInput1, class TOutput=TInput1>
class SubM
{
public:
  typedef typename NumericTraits< TInput1 >::AccumulateType AccumulatorType;
  SubM() {};
  ~SubM() {};
  void SetVal(TInput1 i) { m_Val = i; }

  bool operator!=( const SubM & ) const
  {
    return false;
  }
  bool operator==( const SubM & other ) const
  {
    return !(*this != other);
  }
  inline TOutput operator()( const TInput1 & A, const TInput2 & B)
  {
    if (B) return static_cast<TOutput>(m_Val);
    
    return static_cast<TOutput>( A );
  }
private:
  TInput1 m_Val;
}; 

}
template <class TInputImage1, class TInputImage2, class TOutputImage=TInputImage1>
class ITK_EXPORT SubstituteMaskImageFilter :
    public
BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                         Functor::SubM< 
  typename TInputImage1::PixelType, 
  typename TInputImage2::PixelType,
  typename TOutputImage::PixelType>   >


{
public:
  /** Standard class typedefs. */
  typedef SubstituteMaskImageFilter  Self;
  typedef BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                                   Functor::SubM< 
    typename TInputImage1::PixelType, 
    typename TInputImage2::PixelType,
    typename TOutputImage::PixelType>   
  >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(SubstituteMaskImageFilter, 
               BinaryFunctorImageFilter);

  void SetVal(typename TInputImage1::PixelType val)
  {
    this->GetFunctor().SetVal(val);
    this->Modified();
  }


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
 //  itkConceptMacro(Input1Input2OutputAdditiveOperatorsCheck,
//     (Concept::AdditiveOperators<typename TInputImage1::PixelType,
//                                 typename TInputImage2::PixelType,
//                                 typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  SubstituteMaskImageFilter() {}
  virtual ~SubstituteMaskImageFilter() {}

private:
  SubstituteMaskImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
