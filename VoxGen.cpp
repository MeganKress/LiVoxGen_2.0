/*
VoxGen, 2015

Megan Kress <kressmeg@msu.edu>

This program takes the input file VoxGenPar and its
user-defined categories.

The program may output multiple files:
    - Histogram Flat File
            (see VoxGenFlatFileFormats text file)


*/

#include <stdio.h>
#include <time.h>       /* time */
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstddef>
#include <string>
#include <set>
#include <cmath>
#include <exception>
#include <fstream>  // std::ifstream
#include <liblas/liblas.hpp>
#include "kdtree.h"
#include <iomanip>
#include <map>
#include <math.h>
#include "vgpar.h"
#include "point.h"
#include "voxel.h"
#include "voxcol.h"
#include "voxdata.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>

#include <gdal/cpl_conv.h> // for CPLMalloc()
#include <gdal/gdal_priv.h>
//#include <gdal/gdal_alg.h>
//#include <gdal/gdalwarper.h>
//#include <gdal/ogr_spatialref.h>
//#include <gdal/cpl_conv.h>
//#include <gdal/cpl_minixml.h>

using namespace std;

namespace fs = boost::filesystem;

namespace std
{
    template<> struct less<Point>
    {
       bool operator() (const Point& lhs, const Point& rhs)
       {
           return lhs.x < rhs.x || lhs.y < rhs.y || lhs.z < rhs.z;
       }
    };

}

bool exists(string file)
{
    ifstream infile (file.c_str());
    return infile.good();
};



//Suppress output: https://bbs.archlinux.org/viewtopic.php?id=79378

int main(int argc, char *argv[])
{



    time_t timerI;

    time(&timerI);

    string vgfile;
    bool  fileEntry = false;
    string enteredFile;

    std::streambuf* cout_sbuf;

    bool quiet = false;

    if (argc > 3)
    {
        cout_sbuf = std::cout.rdbuf(); // save original sbuf
        std::ofstream   fout("/dev/null");
        std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'

        fileEntry = true;
        vgfile = argv[1];
        enteredFile = argv[2];
        quiet = true;
    }
    else if (argc > 2)
    {
      fileEntry = true;
      vgfile = argv[1];
      enteredFile = argv[2];
    }
    else if (argc > 1)
    {
        vgfile = argv[1];
    }else {
        cout << "\nEnter /path/to/VoxGenPar: ";
        cin >> vgfile;
        vgfile = vgfile.c_str();
    }

    cout << '\n' << vgfile << '\n';
    vgpar par (vgfile);

    vector<string> filesInDir;

    string fORd = par.getString("fileORdir");

    if(fileEntry)
    {
      filesInDir.push_back(enteredFile);

    }else if(fORd[0] == 'd'){
    //**********************************
    // Boost Filesystem Stuff
    // Adapted from
    // http://www.boost.org/doc/libs/1_31_0/libs/filesystem/example/simple_ls.cpp
    //***********************************

    string inDir = par.getString("inDir");

    fs::path full_path( inDir );

    unsigned long file_count = 0;
    unsigned long dir_count = 0;

    if ( !fs::exists( full_path ) )
  {
    std::cout << "\nDirectory not found";
    return 1;
  }

    if ( fs::is_directory( full_path ) )
  {
    std::cout << "\nIn directory: "
              << full_path << "\n\n";
    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( full_path );
          dir_itr != end_iter;
          ++dir_itr )
    {
      try
      {
        if ( fs::is_directory( *dir_itr ) )
        {
          ++dir_count;
        }
        else
        {
          filesInDir.push_back(dir_itr->path().c_str());
          ++file_count;
        }
      }
      catch ( const std::exception & ex )
      {
        cout << "Error";
      }
    }
    std::cout << "\n" << file_count << " files\n"
              << dir_count << " directories\n";
  }
  else // must be a file
  {
    std::cout << "\nFound File (not directory): " << full_path << "\n";
  }

    //*******************************
    // End Boost Filesystem Stuff
    //*******************************
    }else
    {
        string selectedFile = par.getString("inFile");
        filesInDir.push_back(selectedFile);
    }

    for(unsigned int ii = 0; ii < filesInDir.size(); ++ii) {

    cout << "\n***********************************\n";

    string inFile = filesInDir.at(ii);

    if(par.getString("routput")[0] == 'y') if(!exists(par.getString("outDirectory"))){
        cout << "\nR output directory does not exist. Exiting...\n";
        return 0;
    }
    if(par.getString("hff")[0] == 'y') if(!exists(par.getString("hffDir"))){
        cout << "\nHistogram flat file output directory does not exist. Exiting...\n";
        return 0;
    }
    if(par.getString("lmff")[0] == 'y') if(!exists(par.getString("lmffDir"))){
        cout << "\nLidar metrics flat file output directory does not exist. Exiting...\n";
        return 0;
    }

    ifstream ifs;
    ifs.open( inFile.c_str(), ios::in | ios::binary );

    liblas::Reader reader = liblas::Reader(ifs);

    liblas::Header const& header = reader.GetHeader();

//    cout << '\n' << "Signature: " << header.GetFileSignature() << '\n';
//    cout << "Points count: " << header.GetPointRecordsCount() << '\n';

    liblas::SpatialReference srs = header.GetSRS();

//    cout << "Spatial Reference: " << srs.GetWKT() << '\n';

    cout << "\nMin X: " << setprecision(12) << header.GetMinX() << '\n';
    cout << "Max X: " << setprecision(12) << header.GetMaxX() << '\n';
    cout << "Min Y: " << setprecision(12) << header.GetMinY() << '\n';
    cout << "Max Y: " << setprecision(12) << header.GetMaxY() << '\n';
    cout << "Min Z: " << setprecision(12) << header.GetMinZ() << '\n';
    cout << "Max Z: " << setprecision(12) << header.GetMaxZ() << '\n' << '\n';




    string filter;

    filter = par.getString("filter");

    if(filter[0] == 'y' || filter[0] == 'Y')
    {

    // *************************************
    // Filter by Classification
    // *************************************

    string classChoice;

  std::vector<liblas::Classification> classes;

  classChoice = par.getString("classes");

  for(unsigned int i = 0; i < classChoice.length(); ++i)
    {
        if(classChoice[i] == 'A') break;
        if(classChoice[i] == 'U') classes.push_back(liblas::Classification(1));
        if(classChoice[i] == 'G') classes.push_back(liblas::Classification(2));
        if(classChoice[i] == 'L') classes.push_back(liblas::Classification(3));
        if(classChoice[i] == 'M') classes.push_back(liblas::Classification(4));
        if(classChoice[i] == 'H') classes.push_back(liblas::Classification(5));
        if(classChoice[i] == 'B') classes.push_back(liblas::Classification(6));
        if(classChoice[i] == 'W') classes.push_back(liblas::Classification(9));
    }

    std::vector<liblas::FilterPtr> filters;
    liblas::FilterPtr class_filter = liblas::FilterPtr(new liblas::ClassificationFilter(classes));

    // eInclusion means to keep the classes that match.  eExclusion would
    // throw out those that matched
    class_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(class_filter);

    // ************************************
    // End Filter by Classification
    // ************************************

    // ************************************
    // Filter by return
    // ************************************

    string returnChoice;
    liblas::ReturnFilter::return_list_type returns;


    returnChoice = par.getString("returns");

    for(unsigned int i = 0; i < returnChoice.length(); ++i)
    {
        if(classChoice[i] == 'A') break;
        if(classChoice[i] == '1') returns.push_back(1);
        if(classChoice[i] == '2') returns.push_back(2);
        if(classChoice[i] == '3') returns.push_back(3);
        if(classChoice[i] == '4') returns.push_back(4);
        if(classChoice[i] == '5') returns.push_back(5);
    }


    liblas::FilterPtr return_filter = liblas::FilterPtr(new liblas::ReturnFilter(returns, false));
    return_filter->SetType(liblas::FilterI::eInclusion);
    filters.push_back(return_filter);

    // *************************************
    // End filter by return
    // *************************************

//
//     *************************************
//     Filter by Scan Angle
//     *************************************

    double angleChoice;
    liblas::ContinuousValueFilter<double>::filter_func f;
    liblas::ContinuousValueFilter<double>::compare_func c;
//    cout << "\n\nSelect Maximum Scan Angle.";
//    cout << "\nAngle: ";

    //cin >> angleChoice;

    angleChoice = par.getNum("angle");

    f = &liblas::Point::GetScanAngleRank;
    c = std::less_equal<double>();

    liblas::FilterPtr angle_filterptr = liblas::FilterPtr(new liblas::ContinuousValueFilter<double>(f, angleChoice, c));

    angle_filterptr->SetType(liblas::FilterI::eInclusion);

    filters.push_back(angle_filterptr);

    c = std::greater_equal<double>();
    angleChoice *= -1;

    liblas::FilterPtr angle_filterptr2 = liblas::FilterPtr(new liblas::ContinuousValueFilter<double>(f, angleChoice, c));

    angle_filterptr2->SetType(liblas::FilterI::eInclusion);

    filters.push_back(angle_filterptr2);


    // ***************************************
    // End Filter by Scan Angle
    // ***************************************


   reader.SetFilters(filters);
    }


    kdtree *kd = kd_create(3);
    struct kdres *results;
    cout << "\nAdding points to tree..." << '\n';

   int totalPoints = 0;

   vector<Point> allPoints;



   if(par.getString("hff")[0] == 'y')
    {
      while (reader.ReadNextPoint())
      {
        liblas::Point const& p = reader.GetPoint();

        double point[3];

        point[0] = p.GetX();
        point[1] = p.GetY();
        point[2] = p.GetZ();

//        Point pt (p.GetX(), p.Ge)
//        allPoints.push_back(Point (p.GetX(), p.GetY(), p.GetZ()));

//        cout << p.GetX() << ',' << p.GetY() << ','<< p.GetZ() << '\n';

        ++totalPoints;

        kd_insert3(kd, p.GetX(), p.GetY(), p.GetZ(), NULL );
      };

    }

  if(totalPoints == 0)
  {
    ifs.close();
    if (quiet) cout.rdbuf(cout_sbuf);
    return 0;
  }

  cout << "\nPoints in kdtree: " << totalPoints <<'\n';
  //Set min and max values for each dimension
  double x1 = header.GetMinX();
  double x2 = header.GetMaxX();

  double y1 = header.GetMinY();
  double y2 = header.GetMaxY();

  double z1 = header.GetMinZ();
  double z2 = header.GetMaxZ();

  if(par.getNum("z1") != 0) z1 = par.getNum("z1");

  if(par.getNum("z2") != 0) z2 = par.getNum("z2");

  double zStart;

  double voxelsLo = 0;
  if(par.getNum("zmin") != 0) zStart = par.getNum("zmin");
  else zStart = z1;

  double ht = par.getNum("height");

//  zStart += ht/2.0;
  double z0 = zStart;

  while(zStart <= z1)
  {
    zStart += ht;
    ++voxelsLo;
  }

  zStart -= 2.0*ht;
  --voxelsLo;
  z1 = zStart;

  cout << "\n...Done creating tree.\n";



  double b;
  double h;
  double xcoord;
  double ycoord;
  double xcoord2;
  double ycoord2;


  b = par.getNum("base");


  h = par.getNum("height");


  xcoord = par.getNum("x1");

  if(xcoord == 0) xcoord = x1;

  xcoord2 = par.getNum("x2");

  if(xcoord2 == 0) xcoord2 = x2;

  ycoord = par.getNum("y1");

  if(ycoord == 0) ycoord = y1;


  ycoord2 = par.getNum("y2");

  if(ycoord2 == 0) ycoord2 = y2;


  int numVoxX = ceil((xcoord2 - xcoord) / b);
  ++numVoxX;
  int numVoxY = ceil((ycoord2 - ycoord) / b);
  ++numVoxY;

  int numVoxCols = numVoxX * numVoxY;



  vector<VoxCol> voxColumnsVector;

  int numVox = 0;
  double z = z1;
  while(z <= z2)
  {
    z += ht;
    ++numVox;
  }

  z += ht;
  ++numVox;


  z2 = z;

  b = b/2.0;
  h = h/2.0;


  //Set Radius of Sphere
  double r = sqrt(2*b*b + h*h);

  double x = xcoord;
  double y = ycoord;
  double zcoord = z1;




  int totalnumVox = numVoxCols * numVox;


  int kdPointCount = 0;




  if(par.getString("hff")[0] == 'y')
  {
      cout << "\nNumber of Voxels: " << totalnumVox << '\n';


      cout << "Creating Voxels...\n";

        int hundVox = totalnumVox/100.0;
        int percent = 0;
        int voxCount = 0;

      int voxCt = 0;



      time_t timer;

      time(&timer);

      for(int i = 0; i < numVoxY; ++i)
      {

        for(int j = 0; j < numVoxX; ++j)
        {


          VoxCol voxelColumn;
          voxelColumn.pic = 0;

          for(int k = 0; k < numVox; ++k)
          {


            results = kd_nearest_range3(kd, x, y, zcoord, r);

            int pointsNum = 0;

            Voxel voxel;

            voxel.cX = x;
            voxel.cY = y;
            voxel.cZ = zcoord;
            voxel.d = b;
            voxel.h = h;


            double pos [3] ;

            while( !kd_res_end( results ) ) {

                /* get the position of the current result item */
                kd_res_item( results, pos );

                Point point (pos[0], pos[1], pos[2]);

                if(voxel.inVox(pos[0], pos[1], pos[2])) {
                    voxel.pointsInVox.push_back(point);
                    ++pointsNum;
//                    allPoints.resize(std::remove(allPoints.begin(), allPoints.end(), point) - allPoints.begin());
                }


                /* go to the next entry */
                kd_res_next( results );


          }

            kd_res_free( results );


            voxel.pointNum = voxel.pointsInVox.size();



            voxelColumn.voxels.push_back(voxel);
            voxelColumn.pic += voxel.pointNum;
            kdPointCount += voxel.pointNum;


           zcoord += 2*h;
          }


           voxelColumn.xC = x;
           voxelColumn.yC = y;


           voxelColumn.init(-99);
           if(voxelColumn.pic > 0) voxColumnsVector.push_back(voxelColumn);




          zcoord = z1;
          x = x + 2*b;

        }
          x = xcoord;
          y = y + 2*b;
      }




      cout << "\nDone creating voxels.\n";

       time_t timer2;

      time(&timer2);

        cout << "\nTime Elapsed Creating Voxels: " << timer2 - timer;

        cout << "\nPoints added: " << kdPointCount <<'\n';
        cout << "AllPoints: " << allPoints.size() << '\n';

  } else {

      for(int i = 0; i < numVoxY; ++i)
      {
        if(i == 0) y += b;

        for(int j = 0; j < numVoxX; ++j)
        {

          if(j == 0) x += b;

          VoxCol voxelColumn;

           voxelColumn.xC = x;
           voxelColumn.yC = y;

           voxelColumn.init(-99);
           voxColumnsVector.push_back(voxelColumn);


          x = x + 2*b;

        }
          x = xcoord;
          y = y + 2*b;
      }


    kd_free( kd );

  }




    // ***************************************
    // Parse file name to obtain acquisition
    // and segment identifications.
    // ***************************************



    string fileName = inFile;

    string delimiter = "/";
    string delimiter1 = "_";
    string delimiter2 = ".";
    string delimiter3 = "-";

    string ac, seg, splitNum;


    size_t pos = 0;
    size_t pos1 = 0;
    size_t pos2 = 0;

    while ((pos = fileName.find(delimiter)) != string::npos) {
         fileName.erase(0, pos + delimiter.length());
    }

    ac = par.getString("acquisition");

    if(par.getString("split")[0] == 'y')
    {
        pos = fileName.find(delimiter1);
        seg = fileName.substr(0, pos);
        fileName.erase(0, pos + delimiter1.length());
        pos = fileName.find(delimiter2);
        splitNum = fileName.substr(0, pos);

    } else {

        pos = 0;

        int word = 0;
        while ((pos = fileName.find(delimiter1)) != string::npos) {
            pos = fileName.find(delimiter1);
            if(word != 0) ac += "_";
            ac += fileName.substr(0, pos);
            fileName.erase(0, pos + delimiter1.length());
            ++word;
        }

        pos = fileName.find(delimiter2);
        seg = fileName.substr(0, pos);
        splitNum = "";
    }


    // ***************************************
    // GDAL http://www.gdal.org/gdal_tutorial.html
    //
    // This section is run if the user wishes
    // to create a lidar metrics flat file.
    // ***************************************





    // ***************************************
    // End GDAL Stuff
    // ***************************************



    string outFile;



    // ***************************************
    // Set up end of file names (based on id)
    // ***************************************


    VoxData data;
    data.numCols = numVoxCols;
    data.numHi = numVox;
    data.numPts = kdPointCount;
    data.voxcols = voxColumnsVector;



    data.acqid = ac;

    data.segid = seg;

    if(splitNum.compare("") != 0)
    {
        data.splitid = splitNum;
        splitNum = "-" + splitNum;
    }else {
        data.splitid = "1";
    }


    string id = "_" + ac + "_" + seg + splitNum;


    // ***************************************
    // Histogram flat file creation
    // ***************************************

    double zmax;

    if(par.getNum("zmax") != 0) zmax = par.getNum("zmax");
    else zmax = z2;

    int voxelsHi = 0;
    while(z2 <= zmax)
    {
        z2 += ht;
        ++voxelsHi;
    }

    cout << "Voxels Above: " << voxelsHi << '\n';
    cout << "Voxels Below: " << voxelsLo << '\n';
    cout << "Filled Voxels: " << numVox << '\n';



    string hff = par.getString("hff");
    if(hff[0] == 'y')
    {
        cout << "\n\nWriting data to histogram flat file...\n";
        string hffdir = par.getString("hffDir") + "/" + par.getString("hffTitle");



        //data.toHistFF(hffdir, voxelsHi, voxelsLo - 1, numVox, z0, z2);
        data.toHistFF(hffdir, voxelsHi, voxelsLo - 1, numVox, z0, z2);
        cout << "Done writing data to histogram flat file.\n";
    }


    ifs.close();
    }





  if (quiet) std::cout.rdbuf(cout_sbuf); // restore the original stream buffer

  cout << "Buffer restored\n";

  return(0);
};

