#include "tclap/CmdLine.h"
#include <iomanip>

#include <itkBinaryThresholdImageFilter.h>
#include "ioutils.h"

#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"

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
////////////////////////////////////////////////////////////////////
template <class RImage, class LImage>
typename LImage::Pointer maskJugular(typename RImage::Pointer rawIm,
				     typename LImage::Pointer labelIm)
{
  // use a watershed to segment jugular from everything else
  typename RImage::Pointer gradIm;
  {
  // compute gradient
  typedef typename itk::GradientMagnitudeRecursiveGaussianImageFilter<RImage, RImage> GradMagType;
  typename GradMagType::Pointer gradfilt = GradMagType::New();
  gradfilt->SetInput(rawIm);
  gradfilt->SetSigma(0.5);
  gradIm = gradfilt->GetOutput();
  gradIm->Update();
  gradIm->DisconnectPipeline();
  }

  typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<RImage, LImage> WShedType;
  typename WShedType::Pointer wshed = WShedType::New();
  wshed->SetInput(gradIm);
  wshed->SetMarkerImage(labelIm);
  wshed->SetFullyConnected(false);
  wshed->SetMarkWatershedLine(false);
  // select the vessel and delete the background marker
  typedef typename itk::BinaryThresholdImageFilter<LImage,LImage> ThreshType;
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
  // preprocess the rawIm to remove the jugular
  {
  typename MaskImType::Pointer jugMarkerIm = readImOrient<MaskImType>(CmdLineObj.JugularIm);
  
  typename MaskImType::Pointer Jugular = maskJugular<RawImType, MaskImType>(rawIm, jugMarkerIm);
  writeIm<MaskImType>(Jugular, CmdLineObj.OutputIm);
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

