#include <iostream>
#include <vector>
#include<pcl-1.9/pcl/common/common_headers.h>
#include<pcl-1.9/pcl/io/pcd_io.h>
#include<pcl-1.9/pcl/visualization/pcl_visualizer.h>
#include<pcl-1.9/pcl/visualization/cloud_viewer.h>
#include<pcl-1.9/pcl/console/parse.h>
#include<pcl-1.9/pcl/ModelCoefficients.h>
#include <vtkPlaneSource.h>   
using namespace std;
int PlaneSize = 0;
boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer_;
vtkSmartPointer<vtkPolyData> createPlane(const pcl::ModelCoefficients& coefficients, float scale[2]=nullptr,float center_[3] = nullptr)
{
    vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New ();

    plane->SetNormal (coefficients.values[0], coefficients.values[1], coefficients.values[2]);
    double norm_sqr = coefficients.values[0] * coefficients.values[0]
                      + coefficients.values[1] * coefficients.values[1]
                      + coefficients.values[2] * coefficients.values[2];


    plane->Push(-coefficients.values[3] / sqrt(norm_sqr));
    plane->SetResolution(200, 200);
    plane->Update();

    double pt1[3], pt2[3], orig[3],center[3];
    plane->GetPoint1(pt1);
    plane->GetPoint2(pt2);
    plane->GetOrigin(orig);
    plane->GetCenter(center);

    double _pt1[3], _pt2[3];
    float scale1=1.0;
    float scale2=1.0;
    if( scale!= nullptr )
    {
        scale1 = scale[0];
        scale2 = scale[1];
    }
    if( center_!= nullptr )
    {
        center[0] = center_[0];
        center[1] = center_[1];
        center[2] = center_[2];
    }
    for(int i = 0; i < 3; i++) {
        _pt1[i] = scale1 * (pt1[i] - orig[i]);
        _pt2[i] = scale2 * (pt2[i] - orig[i]);
    }
    for(int i=0; i<3;++i)
    {
        pt1[i] = orig[i] + _pt1[i];
        pt2[i] = orig[i] + _pt2[i];
    }
    plane->SetPoint1(pt1);
    plane->SetPoint2(pt2);
    plane->SetCenter(center);
    //延长origin
    // double _origin[3];
    // for(int i=0; i<3;++i)
    // {
    //    _origin[i] = scale*(orig[i]-pt1[i]);
    // }
    // for(int i=0; i<3;++i)
    // {
    //    orig[i] = pt1[i] + _origin[i];
    // }
    // plane->SetOrigin(orig);

    plane->Update();
    return (plane->GetOutput());
}
class kdtree{
    public:
    kdtree(int axis_,float value_,kdtree* left_,kdtree* right_,vector<int> point_indexs_,bool is_leaf_,float lenth_[3],float center_[3])
    :axis(axis_),value(value_),left(left_),right(right_),point_indexs(point_indexs_),is_leaf(is_leaf_){
        lenth[0] = lenth_[0];
        lenth[1] = lenth_[1];
        lenth[2] = lenth_[2];
        center[0] = center_[0];
        center[1] = center_[1];
        center[2] = center_[2];
    }
    float value;
    int axis;//分割轴;0代表x轴,1代表y轴,2代表z轴
    bool is_leaf;
    vector<int> point_indexs;
    float lenth[3];//长方体
    float center[3];
    kdtree* left;
    kdtree* right;
};
kdtree* build_recursive_kdtree(kdtree* &root,pcl::PointCloud<pcl::PointXYZRGB>::Ptr db,vector<int> point_indexs,int axis,int leaf_size,float lenth_[3] = nullptr,float center_[3] = nullptr){
    if(root == nullptr){
        root = new kdtree(axis,0,nullptr,nullptr,point_indexs,false,lenth_,center_);
    }
    //如果
    if(point_indexs.size() > leaf_size){
        sort(point_indexs.begin(),point_indexs.end(),[&](int a,int b){
            if(axis == 0){
                return db->points[a].x < db->points[b].x;
            }else if(axis == 1){
                return db->points[a].y < db->points[b].y;
            }else if(axis == 2){
                return db->points[a].z < db->points[b].z;
            }
        });
        int middle_index = point_indexs.size()/2 - 1;
        vector<int> PointIndexsLeft(point_indexs.begin(),point_indexs.begin()+middle_index+1);
        vector<int> PointIndexsRight(point_indexs.begin()+middle_index+1,point_indexs.end());
        float leftlenth_[3],rightlenth_[3],leftcenter_[3],rightcenter_[3],center_plane[3],scale[2];
        for(int i = 0;i < 3;i++){
            leftlenth_[i] = lenth_[i];
            rightlenth_[i] = lenth_[i];
            leftcenter_[i] = center_[i];
            rightcenter_[i] = center_[i];
            center_plane[i] = center_[i];
        }
        vector<float> normal{0.,0.,0.};

        switch (axis)
        {
        case 0:
            root->value = (db->points[point_indexs[middle_index]].x + db->points[point_indexs[middle_index+1]].x)/2;//防止切割的超平面穿过点
            leftcenter_[0] = (root->value + center_[0] - 0.5*lenth_[0]) / 2;
            rightcenter_[0] = (root->value + center_[0] + 0.5*lenth_[0]) / 2;
            leftlenth_[0] = root->value - (center_[0] - 0.5*lenth_[0]);
            rightlenth_[0] = (center_[0] + 0.5*lenth_[0] - root->value);
            normal[0] = 1.0;
            center_plane[0] = root->value; 
            scale[0] = lenth_[2];
            scale[1] = lenth_[1];
           // cout << "Plane"<< PlaneSize << ":" << root->value << endl;
            break;
        case 1:
            root->value = (db->points[point_indexs[middle_index]].y + db->points[point_indexs[middle_index+1]].y)/2;
            leftcenter_[1] = (root->value + center_[1] - 0.5*lenth_[1]) / 2;
            rightcenter_[1] = (root->value + center_[1] + 0.5*lenth_[1]) / 2;
            leftlenth_[1] = root->value - (center_[1] - 0.5*lenth_[1]);
            rightlenth_[1] = (center_[1] + 0.5*lenth_[1] - root->value);
            normal[1] = 1.0; 
            center_plane[1] = root->value; 
            scale[0] = lenth_[0];
            scale[1] = lenth_[2];
            // cout << "point1.x" << db->points[point_indexs[middle_index]].x << "point1.y" << db->points[point_indexs[middle_index]].y
            //  << "point2.x" << db->points[point_indexs[middle_index+1]].x << "point2.y" << db->points[point_indexs[middle_index+1]].y << endl;
            // cout << "Plane W:"<< scale[0] << "Plane H:" << scale[1] << endl;
            // cout << "center x:"<< center_plane[0] << "center y:" << center_plane[1] << "center z" << center_plane[2] << endl;
            // cout << "Plane"<< PlaneSize << ":" << root->value << endl;
            break;
        case 2:
            root->value = (db->points[point_indexs[middle_index]].z + db->points[point_indexs[middle_index+1]].z)/2;
            leftcenter_[2] = (root->value + center_[2] - 0.5*lenth_[2]) / 2;
            rightcenter_[2] = (root->value + center_[2] + 0.5*lenth_[2]) / 2;
            leftlenth_[2] = root->value - (center_[2] - 0.5*lenth_[2]);
            rightlenth_[2] = (center_[2] + 0.5*lenth_[2] - root->value);
            normal[2] = 1.0;
            center_plane[2] = root->value;
            scale[0] = lenth_[0];
            scale[1] = lenth_[1]; 
            break;
        }
        pcl::ModelCoefficients coeffs;
        coeffs.values.clear();
        coeffs.values.push_back(normal[0]);
        coeffs.values.push_back(normal[1]);
        coeffs.values.push_back(normal[2]);
        coeffs.values.push_back(0.0);
        auto plane1 = createPlane(coeffs,scale,center_plane);
        viewer_->addModelFromPolyData (plane1, "plane_"+to_string(PlaneSize));
        viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.1, 0.1, 0.1, "plane_"+to_string(PlaneSize), 0);//颜色
        viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY,0.1, "plane_"+to_string(PlaneSize++), 0);//透明度
        
        root->left = build_recursive_kdtree(root->left,db,PointIndexsLeft,(axis+1)%3,leaf_size,leftlenth_,leftcenter_);
        root->right = build_recursive_kdtree(root->right,db,PointIndexsRight,(axis+1)%3,leaf_size,rightlenth_,rightcenter_);
    }else{
            root->is_leaf = true;
    }
    return root;
}
boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud)
{
  // --------------------------------------------
  // -----Open 3D viewer and add point cloud-----
  // --------------------------------------------
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(1, 1, 1);
  pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
  viewer->addPointCloud<pcl::PointXYZRGB>(cloud, rgb, "sample cloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "sample cloud");
  viewer->addCoordinateSystem(1.0);
  viewer->initCameraParameters();
  return (viewer);
}


int main(int argc, char **argv) {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud_ptr (new pcl::PointCloud<pcl::PointXYZRGB>);
    unsigned int PointSize = 100;
    uint8_t r(255), g(15), b(15);
    for (int i = 0;i < PointSize;i++)
    {
        pcl::PointXYZRGB point;
        point.x = 10*(rand()/(RAND_MAX+1.0f));
        point.y = 10*(rand()/(RAND_MAX+1.0f));
        point.z = 10*(rand()/(RAND_MAX+1.0f));
        point.r = (rand()%(255+1));
        point.g = (rand()%(255+1));
        point.b = (rand()%(255+1));
        point_cloud_ptr->points.push_back(point);
       // cout << "point" << to_string(i) << ":" << point.x << "," <<  point.y << "," << point.z << endl;
    }
    point_cloud_ptr->width = (int) point_cloud_ptr->points.size ();
    point_cloud_ptr->height = 1;
    kdtree* root = nullptr;
    vector<int> point_indexs;
    point_indexs.resize(point_cloud_ptr->points.size());
    for(int i = 0 ;i < point_indexs.size();++i){
        point_indexs[i] = i;
    }
   
    float center_[3] = {5,5,5};
    float RectangularSides[3] = {10,10,10};
    viewer_ = rgbVis(point_cloud_ptr);
    build_recursive_kdtree(root,point_cloud_ptr,point_indexs,0,1,RectangularSides,center_);

    while (!viewer_->wasStopped())
    {
        viewer_->spinOnce(100);
       // boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }
    return 0;
}