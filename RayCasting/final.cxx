#define _USE_MATH_DEFINES
#include <math.h>
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
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPNGWriter.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkExtractRectilinearGrid.h>
#include "astro.cxx"
#include <vtkImageData.h>
#include <time.h>



int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}


int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}




int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}



int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}



void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
}



void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
     idx[0] = cellId%(dims[0]-1);
     idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
     idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}
void GetPointVectorValues(int x, int y, int z, const float *F, double *values, const int *dims){
  values[0] = F[ (z*dims[0]*dims[1] + y*dims[0] + x)];
  values[1] = F[ (z*dims[0]*dims[1] + y*dims[0] + (x+1))];
  values[2] = F[ ((z+1)*dims[0]*dims[1] + y*dims[0] + x)];
  values[3] = F[ ((z+1)*dims[0]*dims[1] + y*dims[0] + (x+1))];
  values[4] = F[ (z*dims[0]*dims[1] + (y+1)*dims[0] + x)];
  values[5] = F[ (z*dims[0]*dims[1] + (y+1)*dims[0] + (x+1))];
  values[6] = F[ ((z+1)*dims[0]*dims[1] + (y+1)*dims[0] + x)];
  values[7] = F[ ((z+1)*dims[0]*dims[1] + (y+1)*dims[0] + (x+1))];
  //cout << "values\t" << values[0] <<" " << values[1] << " " << values[2] << endl; 
}


void BoundingBoxForCell(const float *X, const float *Y, const float *Z, const int *dims,
                   int cellId, double *bbox)
{
  if(cellId >= 0 && cellId < GetNumberOfCells(dims)){
                int idx[3];
                GetLogicalCellIndex(idx, cellId, dims);
                bbox[0] = X[idx[0]];
                bbox[1] = X[idx[0]+1];
                bbox[2] = Y[idx[1]];
                bbox[3] = Y[idx[1]+1];
                bbox[4] = Z[idx[2]];
                bbox[5] = Z[idx[2]+1];
                }
        else{
                bbox[0] = -10000;
                bbox[1] = +10000;
                bbox[2] = -10000;
                bbox[3] = +10000;
                bbox[4] = -10000;
                bbox[5] = +10000;
        }
}






double TriLinearInterpolation(const double *pt, const int *dims, const float *X, 
                              const float *Y, const float *Z, const float *F)
{

  double points[8];  
  double val1,val2,val3 = 0;
  int idx[3];
  int i,j,k;
  

  
    for(i=1;i< dims[0];i++){
      if(X[i] > pt[0]){
        idx[0] = i-1;
        
        break;    
      
    }
  }

  
    for(j=1; j< dims[1]; j++){
      if(Y[j] > pt[1]){
        idx[1] = j-1;
        
        break;
      }
    
  }
  
    for(k=1;k < dims[2]; k++){
      if(Z[k] > pt[2]){
        idx[2] = k-1;
        
        break;
      
    }
  }
  if (i >= dims[0] || i ==0 || j >= dims[1] || j ==0){
        return 0;
    }

  else {
    //cout << "ijk\t" << i << " " << j << " " << k << " " << endl;
    int cellId = GetCellIndex(idx,dims);
    double bbox[6];
    BoundingBoxForCell(X,Y,Z,dims,cellId,bbox);
    GetPointVectorValues( idx[0], idx[1], idx[2], F, points, dims);
    double t1 = (pt[0]-bbox[0])/(bbox[1]-bbox[0]);
    double t2 = (pt[1]-bbox[2])/(bbox[3]-bbox[2]);
  //cout << "t1" << t1 << endl;
    double f1 = points[0] + t1*(points[1] - points[0]);
    double f2 = points[4] + t1*(points[5] - points[4]);
    double f3 = points[2] + t1*(points[3] - points[2]);
    double f4 = points[6] + t1*(points[7] - points[6]);
    val1 = f1 + t2*(f2 - f1);
    val2 = f3 + t2 *(f4 - f3);
  //cout << "f1f2 " << f1 <<" "<< f2 << endl;
    
  double t3 = (pt[2]-bbox[4])/(bbox[5]-bbox[4]);
    //cout << "values" << val1 << " " << val2 << endl;
    val3 = val1 + t3*(val2 - val1);
  }
  return val3;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    return image;
}

void get_crossProduct(double vect_A[], double vect_B[], double cross_P[])
 
{
 
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

double get_magnitude(double vec[])
{
  double mag = sqrt (pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));
  return mag;
}

void get_normalVector( double vec[], double mag, double res_vec[]){
    res_vec[0] = vec[0]/mag;
    res_vec[1] = vec[1]/mag;
    res_vec[2] = vec[2]/mag;
}



 
int main( )
{
  int  i,j;                    

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro512.vtk");
    rdr->Update();

    vtkDataSet* dataset = rdr->GetOutput();
    dataset->GetPointData()->SetActiveVectors("grad");
    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);

    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    TransferFunction tf = SetupTransferFunction();
    Camera camera = SetupCamera();
    double look[3];
    double NormLook[3];
    look[0] = (camera.focus[0] - camera.position[0]);
    look[1] = (camera.focus[1] - camera.position[1]);
    look[2] = (camera.focus[2]- camera.position[2]);
    //cout << "look" << look[0] << " " << look[1]<< " " << look[2]<< " " << endl;
    double magn_look = get_magnitude(look);
    get_normalVector(look,magn_look,NormLook);
    //cout << "magn_look" << magn_look << endl;
    int H = 1024;
    int W = 1024;
    vtkImageData *image;
    unsigned char *buffer;
    image = NewImage(H, W);
    buffer = (unsigned char *)image->GetScalarPointer(0,0,0);
    
    float min = 1000.0;
    float maxi =0.0;
    for (i = 0 ; i < 3*H*W ; i++)
        buffer[i] = 0;

    double Cross_r_u[3];
    get_crossProduct(look,camera.up,Cross_r_u);
    double magn_r_u = get_magnitude(Cross_r_u);
    //cout << "ru" << Cross_r_u[0] << " " << Cross_r_u[1] << " " << Cross_r_u[2]<< endl;
    double r_u[3];
    get_normalVector(Cross_r_u,magn_r_u,r_u);

    double Cross_r_v[3];
    get_crossProduct(look,r_u,Cross_r_v);
    //cout << "rv" << Cross_r_v[0] << " " << Cross_r_v[1] << " " << Cross_r_v[2]<< endl;
    double magn_r_v = get_magnitude(Cross_r_v);
    double r_v[3];
    get_normalVector(Cross_r_v,magn_r_v,r_v);
    double r_x[3];
    double r_y[3];
    unsigned char RGB[3];
    double samples = 1000.0;
    double ref_sample = 500.0;
    double del_samp = (camera.far - camera.near)/samples;

    double opacity=0; 
    double total_opacity =0.0;
    bool first_sample = true;
    //cout << "del_sample" << del_samp << endl;
   

    double interm = (2 * tan((camera.angle/2) * M_PI / 180))/W;
    
    r_x[0] = interm * r_u[0];
    r_x[1] = interm * r_u[1];
    r_x[2] = interm * r_u[2];
    //cout << "ray x " << r_x[0] << " "<< r_x[1] << " " << r_x[2] << endl;
    double interm_1 = (2 * tan((camera.angle/2) * M_PI / 180))/H;
    r_y[0] = interm_1 * r_v[0];
    r_y[1] = interm_1 * r_v[1];
    r_y[2] = interm_1 * r_v[2];
    //cout << "ray y " << r_y[0] << " "<< r_y[1] << " " << r_y[2] << endl;
    //cout << "ru" << r_u[0] << " " << r_u[1] << " " << r_u[2] << endl;
    //cout << "rv" << r_v[0] << " " << r_v[1] << " " << r_v[2] << endl;
    double maxOp = 1.0;
    time_t start, end;
    time(&start);
    for(int i =0; i < H;i++){
        for(int j =0;j < W;j++){
            
            
            double r_d[3];
            double temp = ((2.0*i) + 1 -W)/2.0;
            double temp_1 = ((2.0*j) + 1 -H)/2.0;
      //      cout << "temp" << temp << " " << temp_1 << endl;
            r_d[0] = NormLook[0] + (temp * r_x[0]) + (temp_1 * r_y[0]);
            r_d[1] = NormLook[1] + (temp * r_x[1]) + (temp_1 * r_y[1]);
            r_d[2] = NormLook[2] + (temp * r_x[2]) + (temp_1 * r_y[2]);
    //cout << "NormLook " << (look[0]/magn_look) << " " << (look[1]/magn_look) << " " << (look[2]/magn_look)<< endl; 
    //cout << "2" << (temp * r_x[0]) << " "<< temp*r_x[1]  << " " << temp*r_x[2] << endl;
    //cout << "3" << (temp_1 * r_y[0]) << " "<< temp_1 * r_y[1] << " " << temp_1 * r_y[2]<< endl;     
    double magRay = get_magnitude(r_d);
    
    //cout << "r_d " << r_d[0] << " "<< r_d[1] << " " << r_d[2] << endl;
    double ray_pos[3];
   double R,G,B;
    
    first_sample = true;
    R = 0;
    G = 0;
    B = 0;
    float dist = camera.near;
   for(int k = 0;k<samples ;k++){
        
        ray_pos[0] = camera.position[0] + (dist * r_d[0]);
        ray_pos[1] = camera.position[1] + (dist * r_d[1]);
        ray_pos[2] = camera.position[2] + (dist * r_d[2]);
  
        double field = TriLinearInterpolation(ray_pos, dims, X, Y, Z,F);
        tf.ApplyTransferFunction(field, RGB, opacity);
        int bin = tf.GetBin(field);
        
        opacity = 1.0 - pow((1-opacity), ref_sample/samples);
     
        if(first_sample==true){
            total_opacity = opacity;
            R = RGB[0]*opacity;
            G = RGB[1]*opacity;
            B = RGB[2]*opacity;
            first_sample = false;
        }
        else{
            R += (1-total_opacity)*opacity*RGB[0];
            G += (1-total_opacity)*opacity*RGB[1];
            B += (1-total_opacity)*opacity*RGB[2];
            //cout << "RGB\t" << R << "\t" <<  G << "\t" << B << "\t" << endl; 
            total_opacity += (1 - total_opacity)*opacity;
        }
        if(total_opacity >= maxOp)
        {  
        break;
        }
        
  
        dist += del_samp;
       

    }

        
        int offset = 3*(j*W+i);
        buffer[offset + 0] = (unsigned char)R;
        buffer[offset + 1] = (unsigned char)G;
        buffer[offset + 2] = (unsigned char)B;

    }
}
          time(&end);
  double seconds;
  seconds = difftime(end, start);
  cout << "time " << seconds << endl; 
    WriteImage(image, "ray_1024");
}