//***********************************//
// Author: Dr. Reza Darvishi Kamachali.
// All corresponances @Author.
// Only for using after written permission: reza.kamachali@bam.de
//***********************************//

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <random>
#include <cstdlib>
# include <omp.h>


using namespace std;



class Array3D {
    size_t m_width, m_height;
    std::vector<double> m_data;

  public:
    Array3D(size_t x, size_t y, size_t z, double init = 0.):
      m_width(x), m_height(y), m_data(x*y*z, init)
    {}
    double& operator()(size_t x, size_t y, size_t z) {
        return m_data.at(x + y * m_width + z * m_width * m_height);
    }
};




string MakeFileName(string Directory, string NameBase, int TimeStep, string FileExtension)
{
    stringstream converter;
    converter << TimeStep;
    string Count = converter.str();
    int n = Count.size();
    string Number;
    for(int i = 0; i < 8 - n; i++) Number.append("0");
    Number.append(Count);

    return Directory + NameBase + Number + FileExtension;
}




void WriteVTK(Array3D& conc, Array3D& lambda, Array3D& murho, Array3D& muterm23, Array3D& murhocentre, size_t x, size_t y, size_t z, int tStep)
{
    stringstream outbufer;

    int Nx = x-2, Ny = y-2, Nz = z-2;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Composition\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    outbufer << "SCALARS Composition" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";

    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << conc(i, j, k) << "\n";
    }

    outbufer << "SCALARS Density" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";


    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << lambda(i, j, k) << "\n";
    }
    
    outbufer << "SCALARS murho" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";


    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << murho(i, j, k) << "\n";
    }
    
    outbufer << "SCALARS muterm23" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";


    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << muterm23(i, j, k) << "\n";
    }
    
    
    outbufer << "SCALARS murhocentre" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";


    for(int k = 1; k <= Nz; ++k)
    for(int j = 1; j <= Ny; ++j)
    for(int i = 1; i <= Nx; ++i)
    {
        outbufer << murhocentre(0, j, k) << "\n";
    }


    string FileName = MakeFileName("./Output/VTK/", "Composition_", tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}


