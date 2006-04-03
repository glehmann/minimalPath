#include <iomanip>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkTimeProbe.h"
#include "itkLabelOverlayImageFilter.h"

int main(int argc, char * argv[])
{

  const int dim = 2;
  typedef unsigned char PType;
  typedef itk::Image< PType, dim >    IType;
  
  // read the input image
  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  // read the marker image
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();
  
  typedef itk::MinimalPathImageFilter< IType, IType > MinPathType;
  MinPathType::Pointer path = MinPathType::New();
  path->SetInput( reader->GetOutput() );
  path->SetMarkerImage(reader2->GetOutput());
  //itk::SimpleFilterWatcher watcher(path, "path");

  int repeats = atoi(argv[3]);

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
  path->SetLabelChain(order);


  itk::TimeProbe ltime;
  
  for (int i = 0;i<repeats;i++) 
    {
    ltime.Start();
    path->Update();
    ltime.Stop();
    path->Modified();
    }

  std::cout << std::setprecision(3) << "Path time : " 
            << ltime.GetMeanTime() << std::endl;

  const MinPathType::CostVectorType CV = path->GetCosts();
  MinPathType::CostVectorType::const_iterator CVIt;
  for (CVIt = CV.begin();CVIt!=CV.end();CVIt++)
    {
    std::cout << std::setprecision(3) << "Path cost : " << *CVIt << std::endl;
    }

  // set up the overlay here
  typedef itk::RGBPixel<unsigned char> CPType;
  typedef itk::Image< CPType, dim > RGBType;

  typedef itk::LabelOverlayImageFilter<IType, IType, RGBType> OvType;
  OvType::Pointer OvFilt = OvType::New();
  OvFilt->SetInput(reader->GetOutput());
  OvFilt->SetLabelImage(path->GetOutput());
  OvFilt->SetOpacity(1);
  typedef itk::ImageFileWriter< RGBType > RGBWriterType;
  RGBWriterType::Pointer  RGBwriter =  RGBWriterType::New();

  RGBwriter->SetInput(OvFilt->GetOutput());
  RGBwriter->SetFileName(argv[5]);
  RGBwriter->Update();

  // save a marker image
  OvFilt->SetLabelImage(reader2->GetOutput());
  RGBwriter->SetFileName(argv[4]);
  RGBwriter->Update();

  return 0;
}

