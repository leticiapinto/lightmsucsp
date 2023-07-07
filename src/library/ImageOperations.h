Mat create_image(){
    int width = 10;
    int height = 10;
    cout <<"width : "<<width<<", height : "<<height<<endl;
    Mat image = Mat::zeros(height,width, CV_8UC3);

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            Vec3f pixel2; 
            pixel2[0] = 255; //b
            pixel2[1] = 255; //g
            pixel2[2] = 255; //r

            image.at<Vec3b>(i, j) = pixel2;
        }
    }

    Vec3f pixel; 
    pixel[0] = 0; //b
    pixel[1] = 0; //g
    pixel[2] = 0; //r

    for(int i = 2; i < 8; i++ )
    {
        for(int j = 2; j < 8; j++)
        {
            image.at<Vec3b>(i, j) = pixel;
        }
    }

    pixel[0] = 255; //b
    pixel[1] = 255; //g
    pixel[2] = 255; //r

    for(int i = 4; i < 6; i++ )
    {
        for(int j = 4; j < 6; j++)
        {
            image.at<Vec3b>(i, j) = pixel;
        }
    }

    cout<<"cols: "<<image.cols<<", rows : "<<image.rows<<endl;
    return image;
}


Mat create_test_image(){
    int width = 4;
    int height = 4;
    cout <<"width : "<<width<<", height : "<<height<<endl;
    Mat image = Mat::zeros(height,width, CV_8UC3);

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            Vec3f pixel2; 
            pixel2[0] = 255; //b
            pixel2[1] = 255; //g
            pixel2[2] = 255; //r

            image.at<Vec3b>(i, j) = pixel2;
        }
    }

    Vec3f pixel; 
    pixel[0] = 0; //b
    pixel[1] = 0; //g
    pixel[2] = 0; //r

    image.at<Vec3b>(1, 1) = pixel;
    image.at<Vec3b>(1, 2) = pixel;
    image.at<Vec3b>(2, 1) = pixel;
    image.at<Vec3b>(2, 2) = pixel;

    cout<<"cols: "<<image.cols<<", rows : "<<image.rows<<endl;
    return image;
}


Mat generate_image(vector< pair<int, iPair> > Saliency, int width, int height)
{
    Mat image = Mat::zeros(height*2+1,width*2+1, CV_8UC1);
    Mat image_temp = Mat::zeros(height*2+1,width*2+1, CV_8UC1);
    int matrix_image[height*2+2][width*2+2];
    
    for(int i=0; i<height*2+1; i++)
    {
        for(int j=0; j<width*2+1; j++)
        {
            matrix_image[i][j] = 0;
        }
    }
    
    //cout<<"fin matrix with 0"<<endl;

    int max_w =0;
    vector< pair<int, iPair> >::iterator it;    
    int count = 0;
    for (it=Saliency.begin(); it!=Saliency.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        int w = it->first;
        
        if (w!=0)
        {
            //cout<<"w : "<<w;
            //cout<<"width : "<<width<<", height : "<<height<<endl;
            //cout<<", u : "<<u<<", v : "<<v;

            if(w>=max_w)
            {
                max_w =w;
            }
            
            int i_u = int(u/(width));
            int j_u = u%(width);

            //cout<<"       (i_u : "<<i_u<<", j_u : "<<j_u<<" ----\t";
            int i_v = int(v/(width));;
            int j_v =  v%(width); 

            //cout<<"i_v : "<<i_v<<", j_v : "<<j_v<<")  ";
            
            /*
            if(abs(i_u-i_v)==1 || abs(j_u-j_v)==1) //misma fila
            {*/
            int min_i,min_j;
            bool misma_fila =false, misma_columna =false;

            if(j_u==j_v)
            {
                misma_columna =true;
            }

            if(i_u<i_v)
            {
                min_i = i_u;
                min_j = j_u;
                
            }else if(i_u==i_v)
            {
                misma_fila = true;
                min_i = i_u;
                if(j_u<=j_v)
                {
                    min_j = j_u;
                }else{
                    min_j = j_v;
                }
                
            }
            else{
                min_i =i_v;
                min_j =j_v;
            }

            int index_i = 2*min_i +1;
            int index_j = 2*min_j +1;

            //cout<<"\t*** "<<"min_i : "<<min_i<<", min_j : "<<min_j;
            //cout<<", index_i : "<<index_i<<", index_j : "<<index_j;
            if(misma_fila)
            {
                index_j +=1;
                //cout<<" [misma fila] ";
            }
            
            if(misma_columna){
                index_i +=1;
                //cout<<" [misma col] ";
            }

            //cout<<", index_i : "<<index_i<<", index_j : "<<index_j<<endl;
            matrix_image[index_i][index_j] = w;
            /*}*/
        }
    }
 
    for(int i=0; i<height*2+1; i++)
    {
        for(int j=0; j<width*2+1; j++)
        {
            //cout<<matrix_image[i][j]<<"\t";
        }
        //cout<<endl;
    }

    for(int i=0; i<height*2; i++)
    {
        for(int j=0; j<width*2; j++)
        {
            image_temp.at<uchar>(i,j) = int(matrix_image[i][j]);
        }
    }

    //namedWindow( "Display window2", WINDOW_AUTOSIZE );// Create a window for display.
    //imshow( "Display window2", image_temp );  


    for(int i=0; i<height*2+1; i+=2)
    {
        for(int j=0; j<width*2+1; j+=2)
        {
            //recorro los lados para ver que valor tengo como esquina
            int max_corner = matrix_image[i][j];
            if(j-1>=0)
            {
                if(matrix_image[i][j-1]>=max_corner)
                {
                    max_corner = matrix_image[i][j-1];
                }
            }

            if(j+1<width*2+1)
            {
                if(matrix_image[i][j+1]>=max_corner)
                {
                    max_corner = matrix_image[i][j+1];
                }
            }
            if(i-1>=0)
            {
                if(matrix_image[i-1][j]>=max_corner)
                {
                    max_corner = matrix_image[i-1][j];
                }
            }

            if(i+1<=height*2)
            {
                if(matrix_image[i+1][j]>=max_corner)
                {
                    max_corner = matrix_image[i+1][j];
                }
            }

            matrix_image[i][j] = max_corner;
        }
    }
    
    //cout<<endl<<"final matrix"<<endl;
    
    for(int i=0; i<height*2; i++)
    {
        for(int j=0; j<width*2; j++)
        {
            image.at<uchar>(i,j) = int(float(matrix_image[i][j]*255)/float(max_w)); //antigua
            //image.at<uchar>(i,j) = int(matrix_image[i][j]);
            /*
            if(int(float(matrix_image[i][j]*255)/float(max_w))>253)
            {
                image.at<uchar>(i,j) = 255;
            }else
            {
                image.at<uchar>(i,j) = 0;
            }
            */
        }
    }

    return image;
}


Mat show_coloring(vector<int> &vector_color, int sizecolors, int width, int height)
{
    //Mat image = Mat::zeros(height,width, CV_8UC1); gris
    int matrix_image[height][width];

    Mat image = Mat::zeros(height,width, CV_8UC3);
    //cout<<"color size : "<<vector_color.size()<<endl;

    
    vector<Vec3f > colors;
    for(int i = 0; i < sizecolors; i++)
    {
        //cout<<" show_coloring\t"<<i<<endl;
        Vec3f pixel;
        int aux = (rand() % 100) + 10;
        pixel[0] = (rand() % 250) + aux; //b
        pixel[1] = (rand() % 250) + aux; //g +50
        pixel[2] = (rand() % 250) + aux +50; //r
        colors.push_back(pixel);
    }


    //Vec3f pixel2(0, 0 , 0);
    //pixel[0] = 0;
    //pixel[1] = 0;
    //pixel[2] = 0;
    //colors.push_back(pixel2);
    //cout << "width*height+1 : " << (width*height+1) << endl;

    for(int i = 0; i < width*height; i++)
    {
        //cout<<"vector_color["<<i <<"] = "<<vector_color[i]<<endl;
        int temp_i = int(i/width);
        int temp_j = int(i%width);
        //cout<<"temp_i : "<<temp_i<<", temp_j : "<<temp_j<<endl;
        Vec3f pixel; 
        image.at<Vec3b>(temp_i, temp_j) = colors[vector_color[i]];

    }

    //exit(0);
    return image;

}

Mat show_grayscale_coloring(vector<int> &vector_color, int sizecolors, int width, int height)
{
    //Mat image = Mat::zeros(height,width, CV_8UC1); gris
    int matrix_image[height][width];

    //Mat image = Mat::zeros(height,width, CV_8UC3);
    Mat image = Mat::zeros(height,width, CV_8UC1);

    
    //cout<<"color size : "<<vector_color.size()<<endl;

    
    vector<Vec3f > colors;
    for(int i = 0; i < sizecolors; i++)
    {
        //cout<<" show_coloring\t"<<i<<endl;
        Vec3f pixel;
        int aux = (rand() % 100) + 10;
        pixel[0] = (rand() % 250) + aux; //b
        pixel[1] = (rand() % 250) + aux; //g +50
        pixel[2] = (rand() % 250) + aux +50; //r
        colors.push_back(pixel);
    }


    //Vec3f pixel2(0, 0 , 0);
    //pixel[0] = 0;
    //pixel[1] = 0;
    //pixel[2] = 0;
    //colors.push_back(pixel2);
    //cout << "width*height+1 : " << (width*height+1) << endl;

    for(int i = 0; i < width*height; i++)
    {
        //cout<<"vector_color["<<i <<"] = "<<vector_color[i]<<endl;
        int temp_i = int(i/width);
        int temp_j = int(i%width);
        //cout<<"temp_i : "<<temp_i<<", temp_j : "<<temp_j<<endl;
         
        //image[i][j] = vector_color[i]*255;
        //image.at<Vec3b>(temp_i, temp_j) = colors[vector_color[i]];
        image.at<uchar>(temp_i, temp_j) = vector_color[i]*255;

    }

    //exit(0);
    return image;

}
