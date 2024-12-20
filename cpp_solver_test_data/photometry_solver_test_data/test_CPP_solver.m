
%> mex solver_for_matlab/optimize_pose_mex.cpp solver_for_matlab/photometry_minimizer.cpp -I/usr/include/eigen3/

edvo_output_root_path = "/home/chchien/BrownU/research/StereoEdgeOdometry/cpp_solver_test_data/";
reference_img           = importdata(fullfile(edvo_output_root_path, "m_im_reference.txt"));
reference_edge_mask     = importdata(fullfile(edvo_output_root_path, "m_edgeMask.txt"));
current_img             = importdata(fullfile(edvo_output_root_path, "m_im_current.txt"));
current_img_gx          = importdata(fullfile(edvo_output_root_path, "m_gx.txt"));
current_img_gy          = importdata(fullfile(edvo_output_root_path, "m_gy.txt"));
Edges_3D                = importdata(fullfile(edvo_output_root_path, "m_X3D.txt"));

FX = 517.3;
FY = 516.5;
CX = 318.6;
CY = 255.3;
K = [FX, 0, CX; 0, FY, CY; 0, 0, 1];
initial_relative_pose = [0.9996849702688385, -0.01584537505931321, 0.01946495075311933, -0.00448136893902546;
 0.01571022915310382, 0.9998515452497198, 0.007076451208404443, 0.01822712855446191;
 -0.01957419011220217, -0.006768423079097299, 0.9997854967594172, -0.01764391551388047;
 0, 0, 0, 1];
initial_relative_pose = reshape(initial_relative_pose', [16, 1]);

%> Reverse the 2D edge point locations
Edges_2D = Edges_3D ./ Edges_3D(:,3);
Edges_2D = K * Edges_2D';
Edges_2D = Edges_2D';
num_of_edges = size(Edges_2D, 1);

img_width = 640;
img_height = 480;
K = reshape(K', [9,1]);

%> reshape for display
reference_img_ = reshape(reference_img, [640, 480])';
reference_edge_mask_ = reshape(reference_edge_mask, [640, 480])';
current_img_   = reshape(current_img, [640, 480])';

reference_edge_mask = reference_edge_mask ./ 255;

Edges_3D = reshape(Edges_3D', [3*num_of_edges, 1]);

debug = 1;
relative_pose = optimize_pose_mex(current_img, reference_img, reference_edge_mask, Edges_3D, ...
                                  current_img_gx, current_img_gy, K, initial_relative_pose, ...
                                  num_of_edges, img_width, img_height, debug);
relative_pose_to_reference_frame = reshape(relative_pose, [4,4])';

