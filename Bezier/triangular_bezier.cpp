#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "src/bezier.h"

using namespace std;
using namespace Eigen;
void sample_surface(bezier *b,MatrixXd& samples,MatrixXd& Cp,int n){
    samples = MatrixXd::Zero(n,3);
    for (int i=0;i<n;i++){
        double r =  (double)rand()/RAND_MAX;
        double s =  (1-r)*(double)rand()/RAND_MAX;
        double t =  1-r-s;
        samples.row(i) = b -> TriangularBezier(Cp,3,r,s,t);
    }
}

MatrixXi vis_control_polygon(MatrixXd& Cp,int d){
    MatrixXi face(d*d,3);
    int start;
    for (int i=1;i<=d;i++){
        for (int j =0;j<i;j++){
            start = j+((i-1)*i)/2;
            if (j==0){
                
                face.row((i-1)*(i-1)+j)<< start,start+i+1,start+i;
            }
            else{
                face.row((i-1)*(i-1)+2*j-1)<< start,start+i+1,start+i;
                face.row((i-1)*(i-1)+2*j)<<start,start+i,start-1;
            }
        }
    }
    return face;
    
}

int main(int argc, char *argv[])
{

  //MatrixXd control_points;
  MatrixXd sample;
  bezier*  s = new bezier();
  MatrixXd V,CP_matrix,target_PC,pts_color_c,pts_color_t,all_color;
  MatrixXi F,ctrl_poly;
  vector<Vector3d> sample_v;
  s -> read_control("../tri_patch.txt",CP_matrix);
  //control_points = MatrixXd::Random(10,3)/2;
  //sample_surface(s,V,control_points,1000);
  s->gen_tb_surface(CP_matrix,V,F,3,20);
  ctrl_poly = vis_control_polygon(CP_matrix,3);
  MatrixXd all_pts(V.rows()+CP_matrix.rows(),V.cols());
  MatrixXi all_fs(F.rows()+ctrl_poly.rows(),F.cols());
  Eigen::MatrixXd C(all_fs.rows(),3);
  C<<Eigen::RowVector3d(1.0,1.0,0.0).replicate(F.rows(),1),Eigen::RowVector3d(0.2,1.0,0.2).replicate(ctrl_poly.rows(),1);
  all_pts<<V,CP_matrix;
  all_fs<<F,(ctrl_poly.array()+V.rows());
  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);
  RowVector3f last_mouse;
  long sel = -1;
  // Customize the menu
  int PointSize = 10;
  double r_c = 1.0;
  double g_c,b_c = 0.0;
  double r_t = 0.0;
  double g_t,b_t = 0.5;
  string Cs_path,PC_path;  
  // Add content to the default menu window
  Matrix3d color_c,color_t;
  color_c(0,0) = r_c; color_c(1,1) = g_c; color_c(2,2) = b_c;
  color_t(0,0) = r_t; color_t(1,1) = g_t; color_t(2,2) = b_t;



  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Control Points", ImGuiTreeNodeFlags_DefaultOpen))
    {

      ImGui::InputText("CP Path", Cs_path);
      
      // Expose variable directly ...
      ImGui::InputInt("Control point Size", &PointSize, 0, 0);
      ImGui::InputDouble("r_c", &r_c ,0, 0, "%.2f");
      ImGui::InputDouble("g_c", &g_c ,0, 0, "%.2f");
      ImGui::InputDouble("b_c", &b_c ,0, 0, "%.2f");
      //ImGui::InputDouble("r", &r_t ,0, 0, "%.2f");
      //ImGui::InputDouble("g", &g_t ,0, 0, "%.2f");
      //ImGui::InputDouble("b", &b_t ,0, 0, "%.2f");

    
    }

    if (ImGui::CollapsingHeader("Target PC", ImGuiTreeNodeFlags_DefaultOpen)){

      //ImGui::InputText("Target Path", PC_path);
      if (ImGui::Button("Load", ImVec2(-1,0)))
      {
        //s->read_target(PC_path,target_PC);
        s -> read_control(Cs_path, CP_matrix); 
        //s -> convert2vector(fitted_control_points,3,CP_matrix);
        //sample = s -> gen_tp_surface(fitted_control_points,20);
          
        //s -> gen_mesh(V,F,sample);
        
        //plot control points
        /*
        all_pts.resize(CP_matrix.rows()+target_PC.rows(),3);
        all_pts<<CP_matrix,target_PC;
        pts_color_c = MatrixXd::Ones(CP_matrix.rows(),3) * color_c;
        pts_color_t = MatrixXd::Ones(target_PC.rows(),3) * color_t;
        all_color.resize(pts_color_c.rows()+pts_color_t.rows(),3);
        //all_color << pts_color_c,pts_color_t;
        */
        pts_color_c = MatrixXd::Ones(CP_matrix.rows(),3) * color_c;
        //viewer.data().set_mesh(V, F);
        //viewer.data().set_face_based(true);
        viewer.data().point_size = PointSize;
        viewer.data().set_points(CP_matrix,pts_color_c );

      }
      // Add a button
      if (ImGui::Button("Reset color", ImVec2(-1,0)))
      {
      
        
        color_c(0,0) = r_c; color_c(1,1) = g_c; color_c(2,2) = b_c;
        color_t(0,0) = r_t; color_t(1,1) = g_t; color_t(2,2) = b_t;
        pts_color_c = MatrixXd::Ones(CP_matrix.rows(),3) * color_c;
        //pts_color_t = MatrixXd::Ones(target_PC.rows(),3) * color_t;
        //all_color << pts_color_c,pts_color_t;
        viewer.data().point_size = PointSize;
        viewer.data().set_points(CP_matrix, pts_color_c);
      }
    }
  };
  
  const auto & update = [&]()
  {
    //s -> convert2vector(fitted_control_points,3,CP_matrix);
    //sample = s -> gen_tp_surface(fitted_control_points,20);
    //s -> gen_mesh(V,F,sample);
    s->gen_tb_surface(CP_matrix,V,F,3,20);
    all_pts<<V,CP_matrix;
    all_fs<<F,(ctrl_poly.array()+V.rows());
    viewer.data().set_mesh(all_pts, all_fs);
    viewer.data().set_colors(C);
    viewer.data().set_face_based(true);
    viewer.data().set_points(CP_matrix, Eigen::RowVector3d(r_c, g_c, b_c));
    
  };


  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    last_mouse = Eigen::RowVector3f(
      viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
    
    
      // Move closest control point
    Eigen::MatrixXf CP;
    igl::project(
        CP_matrix.cast<float>(),
        viewer.core().view,
        viewer.core().proj, viewer.core().viewport, CP);
    Eigen::VectorXf D = (CP.rowwise()-last_mouse).rowwise().norm();
    sel = (D.minCoeff(&sel) < 30)?sel:-1;
      if(sel != -1)
      {
        last_mouse(2) = CP(sel,2);
        update();
        return true;
    }
    return false;
  };

  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int,int)->bool
  {
    if(sel!=-1)
    {
      Eigen::RowVector3f drag_mouse(
        viewer.current_mouse_x,
        viewer.core().viewport(3) - viewer.current_mouse_y,
        last_mouse(2));
      Eigen::RowVector3f drag_scene,last_scene;
      igl::unproject(
        drag_mouse,
        viewer.core().view,
        viewer.core().proj,
        viewer.core().viewport,
        drag_scene);
      igl::unproject(
        last_mouse,
        viewer.core().view,
        viewer.core().proj,
        viewer.core().viewport,
        last_scene);
      CP_matrix.row(sel) += (drag_scene-last_scene).cast<double>();
      last_mouse = drag_mouse;
      update();
      return true;
    }
    return false;
  };
  
  viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    sel = -1;
    return false;
  };
 
  viewer.data().set_mesh(all_pts, all_fs);
  viewer.data().set_colors(C);

  viewer.data().set_face_based(true);
  //viewer.core().is_animating = true;
  //plot control points
  
  viewer.data().point_size = PointSize;
  viewer.data().set_points(CP_matrix, Eigen::RowVector3d(r_c, g_c, b_c));
  
  viewer.launch();

}
