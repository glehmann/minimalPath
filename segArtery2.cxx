

#include "tclap/CmdLine.h"
#include <iomanip>
// #include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include <itkLabelStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkHistogram.h>


#include "itkExtractImageFilter.h"

#include "ioutils.h"


#include "itkLabelShapeImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkImageDuplicator.h>
#include <itkNeighborhoodAlgorithm.h>

#include <vector>
#include <set>

// stuff for gradient computation
#ifdef MORHPGRAD
#include "itkMorphologicalGradientImageFilter.h"
#else
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkGradientToMagnitudeImageFilter.h>
#endif
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"
#include <itkChangeLabelImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAndImageFilter.h>

#include "itkAbsDiffConstantImageFilter.h"

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
    
    ValueArg<float> radArg("r","radius","radius around the path (mm) ",false, 20, "float");
    cmd.add(radArg);

    UnlabeledMultiArg<int> labvals(std::string("labels"),
 				   std::string("list of label values"),
 				   std::string("integers"));
    cmd.add(labvals);

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.MarkerIm = markArg.getValue();
    CmdLineObj.SubIm = subArg.getValue();
    CmdLineObj.Labels = labvals.getValue();
    
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
/////////////////////////////////////////////////////
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
  
  std::cout << Mn << std::endl;

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

template <class PixType, int dim>
void segArtery(const CmdLineType &CmdLineObj)
{
  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename itk::Image<unsigned char, dim> MaskImType;

  typename RawImType::Pointer rawIm = readImOrient<RawImType>(CmdLineObj.InputIm);
  typename MaskImType::Pointer pointIm = readImOrient<MaskImType>(CmdLineObj.MarkerIm);

  // compute a cost image - carotid should be dark

  typename RawImType::Pointer costIm = computeCostIm<RawImType, MaskImType>(rawIm, pointIm, CmdLineObj.Labels);
  writeIm<RawImType>(costIm, "/tmp/cost.nii.gz");
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

