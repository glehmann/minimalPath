// A test of artery segmentation
// Inputs are the grey scale image (arteries should be dark), and a
// labelled image containing points along the artery.

#include <iomanip>
// #include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include <itkLabelStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkHistogram.h>

#include "itkvanHerkGilWermanBaseImageFilter.h"

#include "itkExtractImageFilter.h"

#include "ioutils.h"

#if 0
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkNeighborhood.h"
#endif

#include "itkLabelShapeImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkImageDuplicator.h>
#include <itkNeighborhoodAlgorithm.h>


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

const int dim = 3;
typedef unsigned short PType;
typedef signed short SPType;
typedef char LPType;
typedef itk::Image< PType, dim >    IType;
typedef itk::Image< SPType, dim >   SType;
typedef itk::Image< LPType, dim >   LType;

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


// a side filling function
template <class ImType, class RawImType>
typename ImType::Pointer FillSides(typename ImType::Pointer InIm, typename RawImType::Pointer RIm, typename ImType::PixelType value)
{
  typedef typename itk::ImageDuplicator<ImType> idupType;
  typename idupType::Pointer dup;
  dup = idupType::New();
  dup->SetInputImage(InIm);
  dup->Update();
  typename ImType::Pointer result = dup->GetOutput();

  typename ImType::RegionType reg = result->GetLargestPossibleRegion();
  typename ImType::RegionType SS;
  typename ImType::RegionType::SizeType sz = reg.GetSize();
  typename ImType::RegionType::IndexType st = reg.GetIndex();
  
  typedef typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<ImType> FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType FaceListType;
  typedef typename FaceCalculatorType::FaceListType::iterator FaceListTypeIt;

  FaceCalculatorType faceCalculator;

  FaceListType faceList;
  FaceListTypeIt fit;
  typename ImType::RegionType::SizeType ksz;
  ksz.Fill(1);
  faceList = faceCalculator(result, result->GetLargestPossibleRegion(), ksz);

  // iterate over the face list to compute the mean border. We'll use
  // the brighter parts to initialize
  fit = faceList.begin();
  ++fit;
  typedef typename itk::ImageRegionIterator<ImType> ItType;
  typedef typename itk::ImageRegionIterator<RawImType> RItType;
  typedef typename std::vector<typename RawImType::PixelType> PixVec;
  PixVec pv;
//   double total=0;
//   long count = 0;
  for (;fit != faceList.end();++fit)
    {
    const typename ImType::RegionType regT = (*fit);
    RItType rit(RIm, (*fit));
    rit.GoToBegin();
    while(!rit.IsAtEnd())
      {
      typename RawImType::PixelType V = rit.Get();
//       total += (double)V;
      pv.push_back(V);
//       ++count;
       ++rit;
      }
    }

  std::sort(pv.begin(), pv.end());
  // iterate over the face list-ignore the first one
  fit = faceList.begin();
  ++fit;
  // would be much better to use a rank filter
  long pos = long(0.15*pv.size());
  std::cout << "Postion = " << pos << " Value = " << (int)pv[pos] << std::endl;

  typename RawImType::PixelType Thresh = typename RawImType::PixelType(pv[pos]);

  for (;fit != faceList.end();++fit)
    {
    const typename ImType::RegionType regT = (*fit);
    typename ImType::RegionType::SizeType szT = regT.GetSize();
//     if (szT[2] > 1)
//       {
//       // this is a vertical side
      ItType it(result, (*fit));
      RItType rit(RIm, (*fit));
      it.GoToBegin();
      rit.GoToBegin();
      while(!it.IsAtEnd())
	{
	typename RawImType::PixelType V = rit.Get();
	if (V > Thresh) 
	  {
	  it.Set(value);
	  }
	++it;
	++rit;
	}
//       }
    }
  return(result);
}

int main(int argc, char * argv[])
{

  
  // read the input image
//   typedef itk::ImageFileReader< IType > ReaderType;
//   ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( argv[1] );
//   reader->Update();


  // read the marker image
//   typedef itk::ImageFileReader< LType > LReaderType;
//   LReaderType::Pointer reader2 = LReaderType::New();
//   reader2->SetFileName( argv[2] );
//   reader2->Update();

  IType::Pointer raw = readIm<IType>(argv[1]);
  LType::Pointer marks = readIm<LType>(argv[2]);

  // use the markers to compute a typical brightness level for the
  // artery
  typedef itk::LabelStatisticsImageFilter<IType, LType> LabStatsType;
  LabStatsType::Pointer labstats = LabStatsType::New();

  labstats->SetInput(raw);
  labstats->SetLabelInput(marks);
  labstats->Update();
  int radius = atoi(argv[3]);

  typedef itk::MinimalPathImageFilter< IType, LType > MinPathType;
  MinPathType::LabelVectorType order;
  if (argc > 6)
    {
    for (int i = 6; i < argc;i++)
      {
      int v = atoi(argv[i]);
      order.push_back(v);
      }
    } 
  else
    {
    std::cout << "Using labels 1 and 2" << std::endl;
    order.push_back(1);
    order.push_back(2);
    }
  
  float Mn = 0.0;
  float mcount = 0;
  // there are room for changes here. Currently an intersection point
  // will get included twice an
  for (int i = 0;i<order.size();i++)
    {
    if (order[i] != 0)
      {
      Mn += labstats->GetMean(order[i]);
      ++mcount;
      //std::cout << (int)order[i] << "\t" << labstats->GetMean(order[i]) << std::endl;
      }
    }
  Mn /= mcount;
  std::cout << "Mean = " << Mn;

  typedef itk::ShiftScaleImageFilter<IType, SType> ShiftScaleType;
  ShiftScaleType::Pointer shift = ShiftScaleType::New();
  shift->SetInput(raw);
  shift->SetShift(-Mn);
  
  typedef itk::AbsImageFilter<SType, IType> AbsType;
  AbsType::Pointer absF = AbsType::New();
  absF->SetInput(shift->GetOutput());

  MinPathType::Pointer path = MinPathType::New();
  path->SetInput( absF->GetOutput() );
  path->SetMarkerImage(marks);
  //itk::SimpleFilterWatcher watcher(path, "path");

  path->SetLabelChain(order);
  //path->SetDirection(2);
  path->SetUseDistWeights(false);
  path->Update();

  //writeIm<LType>(path->GetOutput(), "/tmp/path.nii.gz");
  const MinPathType::CostVectorType CV = path->GetCosts();
  MinPathType::CostVectorType::const_iterator CVIt;
  for (CVIt = CV.begin();CVIt!=CV.end();CVIt++)
    {
    std::cout << std::setprecision(3) << "Path cost : " << *CVIt << std::endl;
    }


  // compute the bounding box of the path, and clip out a region
  // around it
  typedef itk::LabelShapeImageFilter<LType> LabShapeType;
  LabShapeType::Pointer labshape = LabShapeType::New();
  labshape->SetInput(path->GetOutput());
  labshape->Update();
  LabShapeType::BoundingBoxType bb = labshape->GetBoundingBox(order[0]);
//  std::cout << "BB = " << bb << std::endl;
  LabShapeType::BoundingBoxType::SizeType bbs = bb.GetSize();
  LType::SpacingType S = labshape->GetOutput()->GetSpacing();
  LType::SizeType SZ = labshape->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout << "Spacing = " << labshape->GetOutput()->GetSpacing() << std::endl;
  std::cout << "Size = " << SZ << std::endl;
  
  int R = (int)((float)radius/(float)S[0]);
  LabShapeType::BoundingBoxType::IndexType bbi = bb.GetIndex();
  bb.Print(std::cout);

  bbs[0]+= 2*R + 1;
  bbi[0]-= R;
  R = (int)((float)radius/(float)S[1]);
  bbs[1]+= 2*R + 1;
  bbi[1]-= R;

  for (int k = 0;k<2;k++)
    {
    bbi[k] = std::max(0, int(bbi[k]));
    bbs[k] = std::min(int(SZ[k] - bbi[k]), int(bbs[k]));
    }
  bb.SetSize(bbs);
  bb.SetIndex(bbi);

  bb.Print(std::cout);
 
  // now clip out the region we'll operate on with the watershed
  typedef itk::ExtractImageFilter<LType, LType> LabCropType;
  LabCropType::Pointer labcrop = LabCropType::New();
  labcrop->SetExtractionRegion(bb);
  labcrop->SetInput(labshape->GetOutput());
  // clip out ROI from image
  typedef itk::ExtractImageFilter<IType, IType> ImCropType;
  ImCropType::Pointer imcrop = ImCropType::New();
  imcrop->SetExtractionRegion(bb);
  imcrop->SetInput(absF->GetOutput());

  // Fill the sides of the region - these will be markers
  labcrop->Update();
  imcrop->Update();
  LType::Pointer p = FillSides<LType, IType>(labcrop->GetOutput(), 
					     imcrop->GetOutput(), 
					     order[0]+1);
  //writeIm<IType>(imcrop->GetOutput(), "/tmp/weight.nii.gz");
#ifdef MORPHGRAD
  typedef itk::MorphologicalGradientImageFilter<IType, IType> GradFiltType;
  GradFiltType::Pointer grdMag = GradFiltType::New();
  GradFiltType::KernelType kern;
  kern.SetRadius(1);
  kern.CreateKernel();
  grdMag->SetKernel(kern);
  grdMag->SetInput(imcrop->GetOutput());
#else
  // compute a gradient of the control image
  typedef itk::Image< itk::CovariantVector< itk::NumericTraits< PType>::RealType, dim >, dim > GradImType;
  typedef itk::GradientRecursiveGaussianImageFilter<IType,GradImType> GradFiltType;
  typedef itk::GradientToMagnitudeImageFilter<GradImType, IType> GradMagType;
  GradFiltType::Pointer grd = GradFiltType::New();
  GradMagType::Pointer grdMag = GradMagType::New();

  grd->SetInput(imcrop->GetOutput());
  grd->SetSigma(0.5);
  grdMag->SetInput(grd->GetOutput());
#endif
  // Now set up the watershed
  typedef itk::MorphologicalWatershedFromMarkersImageFilter<IType, LType> WShedType;
  WShedType::Pointer wshed = WShedType::New();
  wshed->SetInput(grdMag->GetOutput());
  wshed->SetMarkerImage(p);
  wshed->SetFullyConnected(false);
  wshed->SetMarkWatershedLine(false);

  // select the vessel and delete the background marker
  typedef itk::ChangeLabelImageFilter<LType, LType> SelectType;
  SelectType::Pointer selectArt = SelectType::New();
  selectArt->SetInput(wshed->GetOutput());
  selectArt->SetChange(order[0]+1, 0);

  writeIm<LType>(selectArt->GetOutput(), argv[4]);
  writeIm<LType>(p, "startmark.nii.gz");

  // reset the crop filter to the input
  imcrop->SetInput(raw);
  
  writeIm<IType>(imcrop->GetOutput(), argv[5]);

  return 0;
}

