

#include "tclap/CmdLine.h"
#include <iomanip>
// #include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include <itkLabelStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkHistogram.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSquareImageFilter.h>

#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include "itkExtractImageFilter.h"

#include "ioutils.h"


#include "itkLabelShapeImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkImageDuplicator.h>
#include <itkNeighborhoodAlgorithm.h>

#include <vector>
#include <set>

// stuff for gradient computation
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"
#include <itkChangeLabelImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAndImageFilter.h>

#include "itkAbsDiffConstantImageFilter.h"
#include "morphutils.h"

template <class pixtype>
class MaxFunctor
{
public:
  MaxFunctor(){}
  ~MaxFunctor(){}
  inline pixtype operator()(const pixtype &A, const pixtype &B)
  {
    return std::max(A, B);
  }
};

typedef class CmdLineType
{
public:
  std::string InputIm, MarkerIm, OutputIm, SubIm;
  std::vector<int> Labels;
  float radius;
  bool morphGrad;
};

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("Artery segmenter ", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> markArg("m","marker","marker image",true,"result","string");
    cmd.add( markArg );

    ValueArg<std::string> subArg("s","subset","subset of input image",true,"result","string");
    cmd.add( subArg );

    ValueArg<std::string> outArg("o","output","output image", true,"result","string");
    cmd.add( outArg );
    
    ValueArg<bool> gradArg("l","linear","linear gradient", false, false,"boolean");
    cmd.add( gradArg );

    ValueArg<float> radArg("r","radius","radius around the path (mm) ",false, 20, "float");
    cmd.add(radArg);

    UnlabeledMultiArg<int> labvals(std::string("labels"),
 				   std::string("list of label values"),
				   true,
 				   std::string("integers"));
    cmd.add(labvals);

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.MarkerIm = markArg.getValue();
    CmdLineObj.SubIm = subArg.getValue();
    CmdLineObj.Labels = labvals.getValue();
    CmdLineObj.radius = radArg.getValue();
    CmdLineObj.morphGrad = !gradArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
/////////////////////////////////////////////////////
// this function computes the absolute difference between average
// brightness of regions defined as markers and the rest of the image
template <class RImage, class LImage>
typename RImage::Pointer computeCostIm(typename RImage::Pointer raw,
				       typename LImage::Pointer marker,
				       const std::vector<int> &labels)
{
  typedef itk::LabelStatisticsImageFilter<RImage, LImage> LabStatsType;
  typename LabStatsType::Pointer labstats = LabStatsType::New();

  labstats->SetInput(raw);
  labstats->SetLabelInput(marker);
  labstats->Update();

  float Mn = 0.0;
  // put the labels in a set so that they don't get counted twice
  typedef std::set<int> IntSet;
  IntSet UniqueLabs;
  for (int i = 0;i<labels.size();i++)
    {
    if (labels[i] != 0)
      {
      UniqueLabs.insert(labels[i]);
      }
    }

  for (IntSet::const_iterator i = UniqueLabs.begin();i != UniqueLabs.end();i++)
    {
    Mn += labstats->GetMean(*i);
    }
  Mn /= UniqueLabs.size();
  
  typedef typename itk::AbsDiffConstantImageFilter<RImage, RImage> AbsDiffType;
  typename AbsDiffType::Pointer AbsDiff = AbsDiffType::New();
  AbsDiff->SetInput(raw);
  AbsDiff->SetVal(Mn);
  typename RImage::Pointer result = AbsDiff->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
/////////////////////////////////////////////////////
template <class TImage>
itk::ImageRegion<TImage::ImageDimension> getBB(typename TImage::Pointer mask, 
					       typename TImage::PixelType V,
					       float padding)
{
  typedef typename itk::LabelShapeImageFilter<TImage> LabShapeType;
  typename LabShapeType::Pointer labshape = LabShapeType::New();
  labshape->SetInput(mask);
  labshape->Update();
  typename LabShapeType::BoundingBoxType bb = labshape->GetBoundingBox(V);

  typename LabShapeType::BoundingBoxType::SizeType bbs = bb.GetSize();
  typename TImage::SpacingType S = labshape->GetOutput()->GetSpacing();
  typename TImage::SizeType SZ = labshape->GetOutput()->GetLargestPossibleRegion().GetSize();
//   std::cout << "Spacing = " << labshape->GetOutput()->GetSpacing() << std::endl;
//   std::cout << "Size = " << SZ << std::endl;

  int R = (int)((float)padding/(float)S[0]);
  typename LabShapeType::BoundingBoxType::IndexType bbi = bb.GetIndex();
//  bb.Print(std::cout);

  // these indexes need to match the orientation set in the reader
  bbs[0]+= 2*R + 1;
  bbi[0]-= R;
  R = (int)((float)padding/(float)S[1]);
  bbs[1]+= 2*R + 1;
  bbi[1]-= R;

  for (int k = 0;k<2;k++)
    {
    bbi[k] = std::max(0, int(bbi[k]));
    bbs[k] = std::min(int(SZ[k] - bbi[k]), int(bbs[k]));
    }
  bb.SetSize(bbs);
  bb.SetIndex(bbi);

//  bb.Print(std::cout);

  return(bb);
}
/////////////////////////////////////////////////////
template <class TImage> 
typename TImage::Pointer doCrop(typename TImage::Pointer Im, typename itk::ImageRegion<TImage::ImageDimension> BB)
{
  typedef typename itk::ExtractImageFilter<TImage, TImage> CropType;
  typename CropType::Pointer cropper = CropType::New();
  cropper->SetExtractionRegion(BB);
  cropper->SetInput(Im);
  typename TImage::Pointer result = cropper->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
/////////////////////////////////////////////////////
template<class TImage> 
typename TImage::Pointer createMarker(typename TImage::Pointer path, float radius)
{
  // dilate the path, invert the result and combine
  typename TImage::Pointer dPath = doDilateMM<TImage>(path, radius);

  typedef typename itk::BinaryThresholdImageFilter<TImage,TImage> ThreshType;
  typename ThreshType::Pointer inverter = ThreshType::New();
  inverter->SetInput(dPath);
  inverter->SetLowerThreshold(1);
  inverter->SetUpperThreshold(1);
  inverter->SetInsideValue(0);
  inverter->SetOutsideValue(2);

  typedef typename itk::MaximumImageFilter<TImage, TImage, TImage> MaxType;
  typename MaxType::Pointer maxfilt = MaxType::New();
  maxfilt->SetInput(inverter->GetOutput());
  maxfilt->SetInput2(path);
  typename TImage::Pointer result = maxfilt->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
/////////////////////////////////////////////////////
template  <class RawIm>
typename RawIm::Pointer upsampleIm(typename RawIm::Pointer input, typename RawIm::SpacingType NewSpacing, int interp=0)
{
  const int dim = RawIm::ImageDimension;
  typedef typename RawIm::PixelType PixelType;

  typedef typename itk::ResampleImageFilter<RawIm, RawIm >  ResampleFilterType;
  typedef typename itk::IdentityTransform< double, dim >  TransformType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  input->Update();

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );
  typedef typename itk::LinearInterpolateImageFunction<RawIm, double >  LInterpolatorType;
  typedef typename itk::NearestNeighborInterpolateImageFunction<RawIm, double >  NNInterpolatorType;

  typename ResampleFilterType::InterpolatorPointerType interpolator;
  switch (interp)
    {
    case 0:
      interpolator = NNInterpolatorType::New();
      break;
    case 1:
      interpolator = LInterpolatorType::New();
      break;
    default:
      std::cout << "Unsupported interpolator" << std::endl;
    }

  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );

  const typename RawIm::SpacingType& inputSpacing = input->GetSpacing();
  typename RawIm::SpacingType spacing;
  typename RawIm::SizeType   inputSize = input->GetLargestPossibleRegion().GetSize();
  typename RawIm::SizeType   size;
  typedef typename RawIm::SizeType::SizeValueType SizeValueType;


  for (unsigned i = 0; i < dim; i++)
    {
    //spacing[i] = inputSpacing[i]/factor;
    float factor = inputSpacing[i]/NewSpacing[i];
    size[i] = static_cast< SizeValueType >( inputSize[i] * factor );
    }
//   std::cout << inputSpacing << NewSpacing << std::endl;
//   std::cout << inputSize << size << input->GetOrigin() << std::endl;
  
  resampler->SetSize( size );
  resampler->SetOutputSpacing( NewSpacing );
  resampler->SetOutputOrigin( input->GetOrigin() );
  resampler->SetInput(input);
  typename RawIm::Pointer result = resampler->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
////////////////////////////////////////////////////////////////////

template <class PixType, int dim>
void segArtery(const CmdLineType &CmdLineObj)
{
  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename itk::Image<unsigned char, dim> MaskImType;

  typename RawImType::Pointer rawIm = readImOrient<RawImType>(CmdLineObj.InputIm);
  typename MaskImType::Pointer pointIm = readImOrient<MaskImType>(CmdLineObj.MarkerIm);

  // dodgy use of scoping to make sure that intermediate images get deleted
  typename RawImType::Pointer costIm, gradIm;
  {
  // compute a cost image - carotid should be dark
  typename RawImType::Pointer basecostIm = computeCostIm<RawImType, MaskImType>(rawIm, pointIm, CmdLineObj.Labels);
  // now compute a gradient of the raw image
  if (CmdLineObj.morphGrad)
    {
    std::cout << "Morph gradient" << std::endl;
    gradIm = doGradientOuterMM<RawImType>(rawIm, 0.5);
    }
  else
    {
    typedef typename itk::GradientMagnitudeRecursiveGaussianImageFilter<RawImType, RawImType> GradMagType;
    typename GradMagType::Pointer gradfilt = GradMagType::New();
    gradfilt->SetInput(rawIm);
    gradfilt->SetSigma(0.5);
    gradIm = gradfilt->GetOutput();
    gradIm->Update();
    gradIm->DisconnectPipeline();
    }

  // combine the gradient image and the cost image - use a max
  // function for want of anything better. Might need to be careful
  // about this if gradient changes

  typedef typename itk::SquareImageFilter<RawImType, RawImType> SquareType;
  typename SquareType::Pointer sqrfilt = SquareType::New();
  sqrfilt->SetInput(gradIm);
  
  // set low gradients to zero
  typedef typename itk::ThresholdImageFilter<RawImType> ThreshType;
  typename ThreshType::Pointer thresh = ThreshType::New();

  thresh->SetInput(sqrfilt->GetOutput());
  thresh->SetLower(5); // WARNING - magic number
  thresh->SetUpper(1000);
  thresh->SetOutsideValue((typename RawImType::PixelType)0.0);

  typedef typename itk::AddImageFilter<RawImType, RawImType, RawImType> AddType;
  typename AddType::Pointer addfilt = AddType::New();
  
  addfilt->SetInput(thresh->GetOutput());
  addfilt->SetInput2(basecostIm);
  costIm = addfilt->GetOutput();
  costIm->Update();
  costIm->DisconnectPipeline();
  }
  writeIm<RawImType>(costIm, "/tmp/cst.nii.gz");
  typedef typename itk::MinimalPathImageFilter<RawImType, MaskImType> MinPathType;
  typename MinPathType::Pointer path = MinPathType::New();


  // convert the labels to the correct type
  std::vector<typename MaskImType::PixelType> LVec;
  for (unsigned i = 0; i < CmdLineObj.Labels.size(); i++)
    {
    LVec.push_back(static_cast<typename MaskImType::PixelType>(CmdLineObj.Labels[i]));
    }
  path->SetInput(costIm);
  path->SetMarkerImage(pointIm);
  path->SetLabelChain(LVec);
  path->SetUseDistWeights(false);
  typename MaskImType::Pointer completePath = path->GetOutput();
  completePath->Update();
  completePath->DisconnectPipeline();

  const typename MinPathType::CostVectorType CV = path->GetCosts();
  typename MinPathType::CostVectorType::const_iterator CVIt;
  for (CVIt = CV.begin();CVIt!=CV.end();CVIt++)
    {
    std::cout << std::setprecision(3) << "Path cost : " << *CVIt << std::endl;
    }

  // resample the input images and marker to improve the spacing along
  // the sparse axis - should be able to do this directly on the
  // cropped images and thereby improve performance. Something went
  // wrong and the image ends up blank
  typename RawImType::SpacingType sp = rawIm->GetSpacing();
  if (sp[2] > 0.75) 
    {
    sp[2] = 0.75;
    typename RawImType::Pointer reRaw = upsampleIm<RawImType>(costIm, sp, 1);
    typename RawImType::Pointer reRaw2 = upsampleIm<RawImType>(rawIm, sp, 1);
    typename MaskImType::Pointer reMark = upsampleIm<MaskImType>(completePath, sp, 0);
    //writeIm<RawImType>(reRaw, "/tmp/h.nii.gz");
    costIm = reRaw;
    completePath = reMark;
    rawIm = reRaw2;
    }

  // clip out a box containing the marker we just created
  itk::ImageRegion<dim> BBox = getBB<MaskImType>(completePath, LVec[0], 30.0);
  typename MaskImType::Pointer cPath = doCrop<MaskImType>(completePath, BBox);
  // doesn't really matter whether the original or cost image is used here
  typename RawImType::Pointer cCost = doCrop<RawImType>(costIm, BBox);
  typename RawImType::Pointer cRaw = doCrop<RawImType>(rawIm, BBox);
  
  // now we need to turn the path image into a marker
  typename MaskImType::Pointer finalMarker = createMarker<MaskImType>(cPath, CmdLineObj.radius);
  // now compute a gradient of the raw image - recomputing instead of
  // cropping in case we have resampled
  if (CmdLineObj.morphGrad)
    {
    std::cout << "Morph gradient" << std::endl;
    gradIm = doGradientOuterMM<RawImType>(cRaw, 0.5);
    }
  else
    {
    typedef typename itk::GradientMagnitudeRecursiveGaussianImageFilter<RawImType, RawImType> GradMagType;
    typename GradMagType::Pointer gradfilt = GradMagType::New();
    gradfilt->SetInput(cRaw);
    gradfilt->SetSigma(0.5);
    gradIm = gradfilt->GetOutput();
    gradIm->Update();
    gradIm->DisconnectPipeline();
    }

  // now apply the watershed

  typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, MaskImType> WShedType;
  typename WShedType::Pointer wshed = WShedType::New();
  wshed->SetInput(gradIm);
  wshed->SetMarkerImage(finalMarker);
  wshed->SetFullyConnected(false);
  wshed->SetMarkWatershedLine(false);
  // select the vessel and delete the background marker
  typedef typename itk::BinaryThresholdImageFilter<MaskImType,MaskImType> ThreshType;
  typename ThreshType::Pointer selector = ThreshType::New();
  selector->SetInput(wshed->GetOutput());
  selector->SetLowerThreshold(1);
  selector->SetUpperThreshold(1);
  selector->SetInsideValue(1);
  selector->SetOutsideValue(0);
  writeIm<MaskImType>(selector->GetOutput(), CmdLineObj.OutputIm);
  writeIm<RawImType>(doCrop<RawImType>(rawIm, BBox), CmdLineObj.SubIm);
  writeIm<MaskImType>(finalMarker, "/tmp/marker.nii.gz");
  writeIm<RawImType>(gradIm, "/tmp/grad.nii.gz");
  writeIm<RawImType>(cCost, "/tmp/cost.nii.gz");
}

/////////////////////////////////////////////////////

int main(int argc, char * argv[])
{
  
  const int dim = 3;
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  typedef itk::Image <unsigned char, dim> ImType;
  typedef itk::ImageFileReader<ImType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(CmdLineObj.InputIm.c_str());
//   reader->GetImageIO()->ReadImageInformation();
  reader->Update();
  const std::type_info &pixelDataType =
    reader->GetImageIO()->GetComponentTypeInfo();

  if ( pixelDataType == typeid(unsigned char))
    {
    segArtery<unsigned char, dim>(CmdLineObj);
    }
  else if (pixelDataType == typeid( float ))
    {
    segArtery<float, dim>(CmdLineObj);
    }
  else
    {
    std::cout << "Unsupported type" << std::endl;
    }


  return EXIT_SUCCESS;
}

