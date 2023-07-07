

#include <iostream>

/* 
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <deque>
#include <utility> 
#include <cmath>


#include "library/LElement.h"
#include "library/Graph.h"
#include "library/ImageOperations.h"

#include "HGB_Edward/HgbSegInterval.h"
*/

#include "functions.h"


using namespace cv;
using namespace std;
//#include <experimental/filesystem>
#include <dirent.h>

typedef  pair<int, int> iPair;
//vector< pair<int, iPair> > MST;


int main(int argc, char *argv[])
{    
     //./main ../../test_image_files/ patient001_frame01_2 ../../test_graph_files/patient001_frame01_2.graph ../../leti/patient001_frame01_2/ _4.graph ../../leti/patient001_frame01_2/
    

    /*
    string _patient = argv[1]; // patient001_frame01_2
    int cutleveltree = atoi(argv[2]);
    int typeComputeD = atoi(argv[3]);
    int typeHGB = atoi(argv[4]);

    string _imagedir; //RESOURCES_PATH) + "image_files";
    _imagedir += RESOURCES_PATH;
    _imagedir += "image_files/";
    string _graphdir;
    _graphdir += RESOURCES_PATH;
    //_graphdir += "graph_files/";
    _graphdir += "1_graph_square_files/";
    //string _graphfile2 = "_4.graph";
    string _testdir = RESULTS_PATH;

    string _fileresults = _testdir + _patient + "/";
    string _graphfile = _graphdir + _patient + ".graph";
    cout << endl << _graphfile << endl;
    
    string _resultname = _patient + ".png";
    
    string _imagefile = _imagedir+""+_patient+".png";
    */
    //_imagefile = "../prueba1_0.png";

    //Mat image = imread( _filename.c_str() , IMREAD_COLOR);   // Read the image file  // antes de reducir el command line
    
    
    //Mat image = imread( _imagefile, IMREAD_COLOR);   // Read the image file  // antes de reducir el command line  //borrando temporalmente

    //Mat final_image = image.clone();

    //Mat final_image_qfz = ProcessGraph( 0, image, _graphfile, _fileresults + "_final_leti.graph", 12 );  // codigo propio para procesar grafo  //0 : readGraph, 1: saveGraph
    //Mat final_image_qfz = ProcessGraph( 1, image, _graphfile, _fileresults + "_final_qfz.graph", 12 );  // codigo propio para procesar grafo  //0 : readGraph, 1: saveGraph
    
    /*
    Mat final_image_sm_1 = sm_code(_fileresults+"1"+_graphfile2, 200); // fijo   //sm_code(../../leti/patient001_frame01_2/1_4.graph, 200);
    Mat final_image_sm_3 = sm_code(_fileresults+"3"+_graphfile2, 1300); // fijo
    Mat final_image_sm_5 = sm_code(_fileresults+"5"+_graphfile2, 218); //fijo
    Mat final_image_sm_7 = sm_code(_fileresults+"7"+_graphfile2, 590); // fijo
    */    

    
    //Mat final_imageHGB = 
    //HGB_Edward(_graphfile, _imagefile, typeHGB, _fileresults+ "cut" + to_string(cutleveltree) +  "_ComputeD" + to_string(typeComputeD) + "_hgb" +  to_string(typeHGB) + ".graph", cutleveltree, typeComputeD); //fijo
    
    
    //imwrite( _fileresults + "cut" + to_string(cutleveltree) +  "_ComputeD" + to_string(typeComputeD) + "_hgb" +  to_string(typeHGB) + "_"  +_resultname, final_imageHGB );
    //cout << _fileresults + "cut" + to_string(cutleveltree) +  "_ComputeD" + to_string(typeComputeD) + "_hgb" +  to_string(typeHGB) + "_"  +_resultname << endl;


    //Mat final_image2 = generate_image(final_image.cols,final_image.rows);      
    //Mat cuadrado = create_image();

    //imwrite( _fileresults+ _resultname, final_image );

    //imwrite( _fileresults+ "qfz_" +_resultname, final_image_qfz );
    /*
    cout<<"guardando para sm" <<endl;  

    imwrite( _fileresults+ "sm_1_"   +_resultname, final_image_sm_1 );
    imwrite( _fileresults+ "sm_3_"   +_resultname, final_image_sm_3 );
    imwrite( _fileresults+ "sm_5_"   +_resultname, final_image_sm_5 );
    imwrite( _fileresults+ "sm_7_"   +_resultname, final_image_sm_7 );
    */


    /*
    namedWindow( "Display window2", WINDOW_AUTOSIZE );// Create a window for display.
    imshow( "Display window2", final_image2 );        // Show our image inside it.
    waitKey(0);
    */

    //string filename = "../patient001_frame01_0.txt";
    //string filename = "/data/leticia/training_text/patient024_frame09_3.txt"; 
    //int op = 9;   // possible options are 9 
    //int level_tree = 5; // possible options are 
    //int typeComputeD = 1; // possible options are

    //Mat final_image = generateboundbox( filename, op, level_tree, typeComputeD );
    //imwrite( "/home/leticia/Documents/Windows_UCSP/results/hgb/image_" + to_string(op) + "_" + to_string(level_tree) + "_" + to_string(typeComputeD) + ".jpg", final_image );

    //level_tree=10;
    //Mat final_image = generateboundbox( filename, op, level_tree, typeComputeD );
    //imwrite( "/home/leticia/Documents/Windows_UCSP/results/hgb/image_" + to_string(op) + "_" + to_string(level_tree) + "_" + to_string(typeComputeD) + ".jpg", final_image );
    //Mat final_image = generateboundbox( filename, op, level_tree=15, typeComputeD );
    //imwrite( "/home/leticia/Documents/Windows_UCSP/results/hgb/image_" + to_string(op) + "_" + to_string(level_tree) + "_" + to_string(typeComputeD) + ".jpg", final_image );
    //cout<<"Termine el programa completo"<<endl;


    cout << argv[1] << argv[2] <<endl;
    string _patient = argv[1];
    string _casenumber    = argv[2];
    //string filename = "/data/leticia/training_text/" + _patient;
    //string filename = "/data/leticia/training_text/patient001_frame01_0.txt";
    string filename = _patient;
    string _case    = _casenumber;
    cout << "filename: " <<filename << ", case: "<< _case << endl;
    
    
    int lst_leveltree[5] = { 1, 5, 10, 15, 20};
    int lst_op[4] = {1, 2, 9, 12};
    int lst_D[4] = {1, 2, 3, 4};

    for(int i=0; i < 5; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int k = 0; k < 4; k++)
            {

                //Mat final_image = generateboundbox2( filename, lst_op[j], lst_leveltree[i], lst_D[k]);
                generateboundbox3( filename, lst_op[j], lst_leveltree[i], lst_D[k]);
                //imwrite( "/data/leticia/results/image_" + _case + "_" + to_string(lst_op[j]) + "_" + to_string(lst_leveltree[i]) + "_" + to_string(lst_D[k]) + ".jpg", final_image );
            }
        }
    }

    
   return 0;
}