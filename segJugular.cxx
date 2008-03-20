#include "tclap/CmdLine.h"
#include <iomanip>

#include <itkBinaryThresholdImageFilter.h>
#include "ioutils.h"
#include "morphutils.h"

#include <itkMaximumImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkChangeLabelImageFilter.h>
#include "itkReconstructionByDilationImageFilter.h"
#include <itkMaskImageFilter.h>

typedef class CmdLineType
{
public:
  std::string InputIm, JugularIm, OutputIm, ReorientedIm;
  std::string targetDir;
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

    ValueArg<std::string> jugArg("j","jugular","jugular image",true,"result","string");
    cmd.add( jugArg );

    ValueArg<std::string> outArg("o","output","output image", true,"result","string");
    cmd.add( outArg );

    ValueArg<std::string> reArg("r","reoriented","reoriented image", true,"result","string");
    cmd.add( reArg );

    ValueArg<std::string> debugArg("d","debugdir","output image directory", false,"/tmp/","string");
    cmd.add( debugArg );
    

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.JugularIm = jugArg.getValue();
    CmdLineObj.ReorientedIm = reArg.getValue();
    CmdLineObj.targetDir = debugArg.getValue() + "/";
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}


/////////////////////////////////////////////////////
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


  for (int i = 0; i < dim; i++)
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
template <class LImage>
typename LImage::Pointer makeMarker(typename LImage::Pointer labelIm)
{
  // select the label value 1 (temporary), dilate, invert, combine
  typedef typename itk::BinaryThresholdImageFilter<LImage,LImage> ThreshType;
  typename ThreshType::Pointer selector = ThreshType::New();
  selector->SetInput(labelIm);
  selector->SetLowerThreshold(1);
  selector->SetUpperThreshold(1);
  selector->SetInsideValue(1);
  selector->SetOutsideValue(0);

  typename LImage::Pointer dilated = doDilateMM<LImage>(selector->GetOutput(), 15);
  typename ThreshType::Pointer invertor = ThreshType::New();
  invertor->SetInput(dilated);
  invertor->SetLowerThreshold(1);
  invertor->SetUpperThreshold(1);
  invertor->SetInsideValue(0);
  invertor->SetOutsideValue(2);

  writeIm<LImage>(invertor->GetOutput(), "/tmp/bg.nii.gz");

  typedef typename itk::MaximumImageFilter<LImage, LImage, LImage> MaxType;
  typename MaxType::Pointer maxfilt = MaxType::New();
  maxfilt->SetInput(labelIm);
  maxfilt->SetInput2(invertor->GetOutput());

  typename LImage::Pointer result = maxfilt->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
////////////////////////////////////////////////////////////////////

template <class RImage, class LImage>
typename LImage::Pointer maskJugular(typename RImage::Pointer rawIm,
				     typename LImage::Pointer labelIm)
{
  typedef typename itk::BinaryThresholdImageFilter<LImage,LImage> ThreshType;

  // apply a morphological reconstruction to the raw image
  typedef typename itk::ReconstructionByDilationImageFilter<RImage, RImage> ReconType;
  typedef typename itk::MaskImageFilter<RImage, LImage, RImage> MaskType;
  typename ThreshType::Pointer jugselect = ThreshType::New();
  jugselect->SetInput(labelIm);
  jugselect->SetLowerThreshold(1);
  jugselect->SetUpperThreshold(1);
  jugselect->SetInsideValue(1);
  jugselect->SetOutsideValue(0);
  writeIm<LImage>(jugselect->GetOutput(), "/tmp/jmarker.nii.gz");

  typename MaskType::Pointer masker = MaskType::New();
  masker->SetInput(rawIm);
  masker->SetInput2(jugselect->GetOutput());
  writeIm<RImage>(masker->GetOutput(), "/tmp/rmarker.nii.gz");
  typename ReconType::Pointer recon = ReconType::New();
  typename RImage::Pointer eRaw = doErodeMM<RImage>(rawIm, 1);

  recon->SetMaskImage(eRaw);
  recon->SetMarkerImage(masker->GetOutput());
  writeIm<RImage>(recon->GetOutput(), "/tmp/recon.nii.gz");
  // use a watershed to segment jugular from everything else
  typename RImage::Pointer gradIm;
  {
  // compute gradient
  typedef typename itk::GradientMagnitudeRecursiveGaussianImageFilter<RImage, RImage> GradMagType;
  typename GradMagType::Pointer gradfilt = GradMagType::New();
  gradfilt->SetInput(rawIm);
  gradfilt->SetSigma(0.75);
  gradIm = gradfilt->GetOutput();
  gradIm->Update();
  gradIm->DisconnectPipeline();
  }

  typename LImage::Pointer marker = makeMarker<LImage>(labelIm);
  writeIm<LImage>(marker, "/tmp/marker.nii.gz");
  writeIm<RImage>(gradIm, "/tmp/gradient.nii.gz");
  typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<RImage, LImage> WShedType;
  typename WShedType::Pointer wshed = WShedType::New();
  wshed->SetInput(gradIm);
  wshed->SetMarkerImage(marker);
  wshed->SetFullyConnected(true);
  wshed->SetMarkWatershedLine(true);
  // select the vessel and delete the background marker
  typename ThreshType::Pointer selector = ThreshType::New();
  selector->SetInput(wshed->GetOutput());
  selector->SetLowerThreshold(1);
  selector->SetUpperThreshold(1);
  selector->SetInsideValue(1);
  selector->SetOutsideValue(0);

  typename LImage::Pointer result = selector->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
////////////////////////////////////////////////////////////////////
template <class PixType, int dim>
void segJugular(const CmdLineType &CmdLineObj)
{
  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename itk::Image<unsigned char, dim> MaskImType;

  typename RawImType::Pointer rawIm = readImOrient<RawImType>(CmdLineObj.InputIm);
  typename MaskImType::Pointer jugMarkerIm = readImOrient<MaskImType>(CmdLineObj.JugularIm);
  typename RawImType::SpacingType sp = rawIm->GetSpacing();

//   if (sp[2] > 0.75) 
//     {
//     sp[2] = 0.75;
//     rawIm = upsampleIm<RawImType>(rawIm, sp, 1);
//     jugMarkerIm = upsampleIm<MaskImType>(jugMarkerIm, sp, 0);
//     }
  // preprocess the rawIm to remove the jugular
  {
  
  
  typename MaskImType::Pointer Jugular = maskJugular<RawImType, MaskImType>(rawIm, jugMarkerIm);
  writeIm<MaskImType>(Jugular, CmdLineObj.OutputIm);
  //writeIm<MaskImType>(jugMarkerIm, "/tmp/marker.nii.gz");  
  }
    
  writeIm<RawImType>(rawIm, CmdLineObj.ReorientedIm);
  return;
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

  if (pixelDataType == typeid( float ))
    {
    segJugular<float, dim>(CmdLineObj);
    }
  else
    {
    std::cout << "Unsupported type" << std::endl;
    }


  return EXIT_SUCCESS;
}

