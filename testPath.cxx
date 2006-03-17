#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMinimalPathImageFilter.h"
#include "itkSimpleFilterWatcher.h"

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
  
  int height = atoi(argv[4]);
  int F = atoi(argv[5]);
  typedef itk::MinimalPathImageFilter< IType, IType > MinPathType;
  MinPathType::Pointer path = MinPathType::New();
  path->SetInput( reader->GetOutput() );
  
  itk::SimpleFilterWatcher watcher(path, "path");
 

  path->Update();

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  std::cout << "writing " << argv[2] << std::endl;
  writer->SetFileName( argv[2] );
  writer->SetInput(path->GetOutput());
  writer->Update();

  

  return 0;
}

