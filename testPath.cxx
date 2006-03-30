#include <iomanip>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkTimeProbe.h"

int main(int, char * argv[])
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

  MinPathType::LabelVectorType order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(3);
  order.push_back(1);
  path->SetLabelChain(order);

//  path->SetStartLabel(1);
//  path->SetEndLabel(2);

  itk::TimeProbe ltime;
  
  for (int i = 0;i<100;i++) 
    {
    ltime.Start();
    path->Update();
    ltime.Stop();
    path->Modified();
    }

  std::cout << std::setprecision(3) << "Path time : " 
            << ltime.GetMeanTime() << std::endl;

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  std::cout << "writing " << argv[3] << std::endl;
  writer->SetFileName( argv[3] );
  writer->SetInput(path->GetOutput());
  writer->Update();

  

  return 0;
}

