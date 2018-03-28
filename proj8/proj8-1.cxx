#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <string>
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include <vtkArrowSource.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkAppendPolyData.h>
#include <vtkAutoInit.h>
#include <vtkLookupTable.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkStreamTracer.h>
#include <vtkPlaneSource.h>
#include <vtkPlane.h>
#include <vtkReverseSense.h>
#include <vtkCutter.h>
#include <vtkGenericClip.h>
#include <vtkSmartPointer.h>
#include "vtkHedgeHog.h"
#include <vtkFieldData.h>
#include <vtkStructuredGridReader.h>
#include <vtkCellData.h>
#include <vtkGlyph3D.h>
#include <vtkLineSource.h>
#include <vtkRungeKutta4.h>
#include <vtkExtractRectilinearGrid.h>
 
int main( )
{
  int  i,j;                    

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");
    rdr->Update();

    vtkDataSet* dataset = rdr->GetOutput();
    dataset->GetPointData()->SetActiveVectors("grad");
    int dims[3];
    vtkSmartPointer<vtkRenderWindow> renderWindow = 
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(800,800);
 
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
 
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // Define viewport ranges
  double xmins[4] = {0,.5,0,.5};
  double xmaxs[4] = {0.5,1,0.5,1};
  double ymins[4] = {0,0,.5,.5};
  double ymaxs[4]= {0.5,0.5,1,1};
  for(unsigned i = 0; i < 4; i++)
    {
    vtkSmartPointer<vtkRenderer> renderer = 
      vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkSphereSource> sphereSource;
    vtkSmartPointer<vtkImageData> volume;
    vtkSmartPointer<vtkContourFilter> cf;
    vtkSmartPointer<vtkContourFilter> cf1;
     vtkSmartPointer<vtkHedgeHog> hedgehog;
    vtkSmartPointer<vtkReverseSense> reverse;
    vtkSmartPointer<vtkReverseSense> reverse1;
    vtkSmartPointer<vtkStreamTracer> streamLine;
    vtkSmartPointer<vtkPlane> plane;
    vtkSmartPointer<vtkPlane> plane1;
    vtkSmartPointer<vtkPlane> plane2;
     vtkSmartPointer<vtkArrowSource> arrow;
    vtkSmartPointer<vtkGlyph3D> glyph3D;
    vtkSmartPointer<vtkLineSource> lineSource;
    vtkSmartPointer<vtkRungeKutta4> RK4;
    vtkSmartPointer<vtkLookupTable> lookupTable;
    vtkSmartPointer<vtkCutter> cutter; 
     vtkSmartPointer<vtkCutter> cutter1;
      vtkSmartPointer<vtkCutter> cutter2; 
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkPolyDataMapper> mapper1;
    vtkSmartPointer<vtkActor> actor1;
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkAppendPolyData> appendFilter;
    renderer->SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i]);
    vtkExtractRectilinearGrid* grid = vtkExtractRectilinearGrid::New();
    double *s;
    lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetTableRange(0.0, 1.0);
    lookupTable->SetHueRange(0.667,0);

    switch(i)
    {
      case 0:
      cf = vtkSmartPointer<vtkContourFilter>::New();
      cf->SetInputConnection(rdr->GetOutputPort());
      cf->SetValue(0,2.5);
      cf1 = vtkSmartPointer<vtkContourFilter>::New();
      cf1->SetInputConnection(rdr->GetOutputPort());
      cf1->SetValue(0,5.0);
      reverse = vtkSmartPointer<vtkReverseSense>::New();
      reverse->SetInputConnection(cf->GetOutputPort());
      reverse->ReverseCellsOn();
      reverse->ReverseNormalsOn();
      reverse1 = vtkSmartPointer<vtkReverseSense>::New();
      reverse1->SetInputConnection(cf1->GetOutputPort());
      reverse1->ReverseCellsOn();
      reverse1->ReverseNormalsOn();
      mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(reverse->GetOutputPort());
      mapper->SetLookupTable(lookupTable);
      mapper->SetScalarRange(0,5);
      mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper1-> SetInputConnection(reverse1->GetOutputPort());
      actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0.0,0.0,0.0);
      renderer->AddActor(actor);
      renderer->ResetCamera();
      actor1 = vtkSmartPointer<vtkActor>::New();
      actor1->SetMapper(mapper1);
      actor1->GetProperty()->SetColor(0.0,0.0,0.0);
      renderer->AddActor(actor1);
      renderer->ResetCamera();

      break;
      case 1:
      dataset->GetPointData()->SetActiveVectors("grad");
      grid->AddInputData(dataset);
      grid->SetSampleRate(5,5,5); 
      grid->Update();
      glyph3D= vtkSmartPointer<vtkGlyph3D>::New();
      arrow = vtkSmartPointer<vtkArrowSource>::New();
      glyph3D->SetSourceConnection(arrow->GetOutputPort());
      glyph3D->SetVectorModeToUseNormal();
      glyph3D->SetInputConnection(grid->GetOutputPort());
      glyph3D->SetVectorModeToUseVector();
      glyph3D->Update();
      mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(glyph3D ->GetOutputPort());
      mapper->SetLookupTable(lookupTable);
      mapper->SetScalarRange(0,6);
      actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      actor->GetProperty()->SetColor(0,0,0);
      renderer->AddActor(actor);
      renderer->ResetCamera(); 

     break;
      case 3: 
      lineSource = vtkSmartPointer<vtkLineSource>::New();
      lineSource->SetPoint1(9.0, 0, 0);
      lineSource->SetPoint2(-9.0, 0, 0);
      lineSource->SetResolution(18);
      lineSource->Update();
      streamLine = vtkSmartPointer<vtkStreamTracer>::New();
      streamLine->SetInputData(dataset);
      streamLine->SetSourceData(lineSource->GetOutput());
      streamLine->SetMaximumPropagation(3000);
      streamLine->SetIntegratorTypeToRungeKutta4();
      mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(streamLine->GetOutputPort());
      s = dataset -> GetScalarRange();
      mapper->SetScalarRange(s);
      actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      renderer->AddActor(actor);
      renderer->ResetCamera();
      renderWindow->AddRenderer(renderer);
      renderWindowInteractor->SetRenderWindow(renderWindow);
      break;
      case 2: 
      lookupTable = vtkSmartPointer<vtkLookupTable>::New();
      lookupTable->SetTableRange(0.0, 1.0);
      lookupTable->SetHueRange(0.667,0);
      lookupTable->Build();
      s = dataset -> GetScalarRange();
      plane = vtkSmartPointer<vtkPlane>::New();
      plane->SetOrigin(0,0,0);
      plane->SetNormal(1,0,0);

      plane1 = vtkSmartPointer<vtkPlane>::New();
      plane1->SetOrigin(0,0,0);
      plane1->SetNormal(0,1,0);

      plane2 = vtkSmartPointer<vtkPlane>::New();
      plane2->SetOrigin(0,0,0);
      plane2->SetNormal(0,0,1);
      cutter = vtkSmartPointer<vtkCutter>::New();
      cutter->SetCutFunction(plane);
      cutter->SetInputData(dataset);
      cutter->Update();
      cutter1 = vtkSmartPointer<vtkCutter>::New();
      cutter1->SetCutFunction(plane1);
      cutter1->SetInputData(dataset);
      cutter1->Update();
      cutter2 = vtkSmartPointer<vtkCutter>::New();
      cutter2->SetCutFunction(plane2);
      cutter2->SetInputData(dataset);
      cutter2->Update();
      appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
      appendFilter->AddInputConnection(cutter->GetOutputPort());
      appendFilter->AddInputConnection(cutter1->GetOutputPort());
      appendFilter->AddInputConnection(cutter2->GetOutputPort());
      appendFilter->Update();
      mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(appendFilter->GetOutputPort());
      mapper->SetLookupTable(lookupTable);

      mapper->SetScalarRange(s);
      actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);
      renderer->AddActor(actor);
      renderer->ResetCamera();
     
  
        break;
        
      default:
        break;
    }
 

    renderWindow->Render();
    renderWindow->SetWindowName("Multiple ViewPorts");
    }
 
  renderWindowInteractor->Start();
 
  
}