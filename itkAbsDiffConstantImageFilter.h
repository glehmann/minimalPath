#ifndef __itkAbsDiffConstantImageFilter_h
#define __itkAbsDiffConstantImageFilter_h

namespace itk
{

/** \class AbsDiffConstantImageFilter
 * \brief Computes the absolute difference between an image and a
 * constant. Can be done with ShiftScale and AbsIamgeFilters.
 */


namespace Function {

template< class TInput, class TOutput>
class AbsDiffConst
{
public:
  AbsDiffConst() {m_Val = 0.0;}
  void SetVal(double i) { m_Val = i; }

  ~AbsDiffConst() {}
  bool operator!=( const AbsDiffConst & ) const
  {
    return false;
  }
  bool operator==( const AbsDiffConst & other ) const
  {
    return !(*this != other);
  }
  inline TOutput operator()( const TInput & A )
  { 
    double diff = (double)A - m_Val;
    const double absdiff = ( diff > 0.0 ) ? diff : -diff;
    return static_cast<TOutput>(absdiff);
  }
private:
  double m_Val;
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT AbsDiffConstantImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage,
                        Function::AbsDiffConst<
  typename TInputImage::PixelType,
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef AbsDiffConstantImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, Function::AbsDiffConst< typename TInputImage::PixelType,  typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  void SetVal(double val)
  {
    this->GetFunctor().SetVal(val);
    this->Modified();
  }


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToDoubleCheck,
    (Concept::Convertible<typename TInputImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToOutputCheck,
    (Concept::Convertible<double, typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  AbsDiffConstantImageFilter() {}
  virtual ~AbsDiffConstantImageFilter() {}

private:
  AbsDiffConstantImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif

