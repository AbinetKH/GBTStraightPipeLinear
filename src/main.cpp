#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <omp.h>
#include <string>
#include <chrono>
#include <iomanip>      // std::setprecision



#include "mesh.h"
#include "LM_connectivity_matrix.h"
#include "ik_tensor.h"
#include "element_stiffness_matrix.h"
#include "boundary_conditions.h"
#include "loading.h"
#include "amplitude_function.h"
#include "displacements.h"

#include <matplot/matplot.h>
#include "gnuplot-iostream.h"


int main()
{
	std::string again="yes";
	do {
	/**
	* Preprocessing
	* Input parameters
	* r(pipe radius), length(total length of the pipe), t(wall thickness of the pipe), E(modules of elasticity), K(plate bending stiffness),Q(plate extensional stiffness),q(area load)
	* G(shear modules), no_beam_elem(number of beam elements),l(element length)
	* cross_sect_ref( number of discretization in the cross-section, only for plotting purpose),
	* Boundary conditions (torsion, axisymmetric, extension, bending, local dofs) at the two end of the pipe 1= fixed and 0=free
	* Initialize vector of GBT mode list
	*/
	int cross_sect_ref = 81, no_beam_elem = 25;
	double r = 500, length = 1000, pi = 3.14159, t = 10, E = 205000, mu = 0.3, q = 1, plot_scale = 50;


	std::cout << "FEM_GBT_linear is a C++ code developed for stress and deformation analysis of straight thin-walled circular pipes based on the Generalized Beam Theory(GBT). \n" << std::endl;
	
	std::cout << "Specify the following parameters or enter to use the values: \n" << std::endl;
    std::cout << "Use consistent unit: ton, mm, s, N, MPa \n" << std::endl;
	std::cout << " Number of cross-sectional refinement [default = 81], \n Number of beam elements [default = 25], \n Pipe length [default = 1000], \n Pipe radius [default = 500], " <<
		" \n Pipe thickness [default = 10], \n Elastic modulus [default = 205000], \n Poisson's ratio [default = 0.3], \n Projected area load  [default = 1 N/mm2], \n Plot scale  [default = 50], \n";
	if (std::cin.peek() == '\n') { 

	}
	else if (!(std::cin >> cross_sect_ref >> no_beam_elem >> length >> r >> t >> E >> mu >> q >> plot_scale)) {
		while (!(std::cin >> cross_sect_ref >> no_beam_elem >> length >> r >> t >> E >> mu >> q >> plot_scale)) {
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
			std::cout << "Oops, the input parameter is invalid. Please try again from the beginning.\n" << std::endl;
		}
	}

	double K = (E * pow(t, 3)) / (12 * (1 - pow(mu, 2)));
	double Q = (E * t) / (1 - pow(mu, 2));
	double G = E / (2 * (1 + mu));
	double l = length / no_beam_elem;
	std::vector<int> torsion_dof = { 1,0 }, axisymmetric_dof = { 1,0 }, extension_dof = { 1,0 }, bending_dof = { 1,0 }, allocal_dof = { 1,0 };
	/*std::vector<std::string> mode_list = { "t|t" , "a|a", "1|1", "2|c", "2|v", "2|u", "3|c", "3|v", "3|u", "4|c", "4|v", "4|u", "5|c", "5|v", "5|u" ,
		"6|c", "6|v", "6|u", "7|c", "7|v", "7|u", "8|c", "8|v", "8|u", "9|c", "9|v", "9|u"};, "10|c", "10|v", "10|u", "11|c", "11|v", "11|u",
		"12|c", "12|v", "12|u","13|c", "13|v", "13|u","14|c", "14|v", "14|u","15|c", "15|v", "15|u","16|c", "16|v", "16|u","17|c", "17|v", "17|u",
	 };*/
	std::vector<std::string> mode_list = { "a|a","1|1", "3|c", "3|v", "3|u", "5|c", "5|v", "5|u", "7|c", "7|v", "7|u" , "9|c", "9|v", "9|u" , "10|c", "10|v", "10|u"};/*, "11|c", "11|v", "11|u",
	"12|c", "12|v", "12|u","13|c", "13|v", "13|u","14|c", "14|v", "14|u","15|c", "15|v", "15|u","16|c", "16|v", "16|u","17|c", "17|v", "17|u",
 };*/
	int length_mode_list = std::size(mode_list);
	std::cout << "mode_list  " << std::size(mode_list) << std::endl;



	/*
	* for Qt input interface to be compeleted ...
	std::string mode_type; // shear_v, shear_u, classical, all
	std::string mode_even_odd; //odd, even, all
	int number_of_modes = 50;
	for (int i = 0; i < number_of_modes; i++) {

	}
	*/

	auto start = std::chrono::high_resolution_clock::now();
	/**
	* mesh generator
	*/
	Eigen::MatrixXd y, z, x;
	RectangularMesh obMesh;
	obMesh.meshStraightPipe(cross_sect_ref, no_beam_elem, r, length, &y, &z, &x);

	//std::cout << x << std::endl;
	//std::cout << x.size() << '/n';
	//std::cout << y << std::endl;
	//std::cout << z << std::endl;

	/**
	* initialize LM connectivity matrix for assembly of the global stiffness matrix
	*/
	LmConnectivityMatrix obLM;
	LmConnectivityMatrix obLMelement;
	Eigen::VectorXd LM_element = obLMelement.openEndPipe(mode_list, 1);
	Eigen::MatrixXd LM_matrix = obLM.openEndPipe(mode_list, no_beam_elem);
	//std::cout << "LM_matrix\n  " << LM_matrix << std::endl;


	/**
	* initialize GBT mode second-order coupling tensor
	*/
	couplingTensors obTensor;
	Eigen::MatrixXd ik_tensor = obTensor.ikTensor(mode_list, LM_element, length_mode_list, mu, t, r, Q, K, G);
	//std::cout << "GBT mode second-order coupling tensor \n  " << ik_tensor << std::endl;


	/**
	* build the element stiffness matrix
	*/
	Eigen::MatrixXd element_matrix;
	GbtElementMatrix obStiffness;
	element_matrix = obStiffness.ikLinear(ik_tensor, l, LM_element.rows());
	//std::cout << "element matrix\n   " << element_matrix << std::endl;

	/**
	* assembel the global stiffness matrix
	*/
	int total_dof = LM_matrix(LM_matrix.rows() - 1, LM_matrix.cols() - 1) + 1;
	std::cout << "total dof   " << total_dof << std::endl;
	Eigen::MatrixXd global_system_stiffness_matrix = Eigen::MatrixXd::Zero(total_dof, total_dof);


	//auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
	for (int i = 0; i < no_beam_elem; i++) {
		global_system_stiffness_matrix(LM_matrix(Eigen::all, i), LM_matrix(Eigen::all, i)) += element_matrix;
		//printf("Thread %d, iteration %d \n", omp_get_thread_num(), i);
	}
	//std::cout << "global system matrix\n  " << global_system_stiffness_matrix << std::endl;

	/**
	* apply the boundary conditions
	*/
	BoundaryCondition obBC;
	Eigen::VectorXd fixedDof = obBC.dofsFixed(mode_list, length_mode_list, no_beam_elem, torsion_dof, axisymmetric_dof, extension_dof, bending_dof, allocal_dof, LM_matrix);
	double penality_vlaue = pow(element_matrix.maxCoeff(), 5);
#pragma omp parallel for 
	for (int i = 0; i < fixedDof.rows(); i++) {
		global_system_stiffness_matrix((int)fixedDof(i), (int)fixedDof(i)) = penality_vlaue;
	}
	//std::cout << "list of fixed dof\n  " << fixedDof << std::endl;
	//std::cout << "global system matrix\n  " << global_system_stiffness_matrix << std::endl;

	/**
	* determine the external load conditions
	* this code is only for projected area load
	* to do: add general type of loading
	*/
	GbtLoading obLoading;
	Eigen::VectorXd element_loading = obLoading.loadVector(mode_list, r, l, q, LM_element, 3 * pi / 2, pi / 2);
	Eigen::VectorXd global_system_load_vector = Eigen::VectorXd::Zero(total_dof);
	//std::cout << "element load vector\n  " << element_loading << std::endl;
#pragma omp parallel for 
	for (int i = 0; i < no_beam_elem; i++) {
		global_system_load_vector(LM_matrix(Eigen::all, i)) += element_loading;
		//printf("Thread %d, iteration %d \n", omp_get_thread_num(), i);
	}
	/**
	* Solver
	* solve for the displacements
	*/

	Eigen::VectorXd global_displ(total_dof);

	global_displ = global_system_stiffness_matrix.llt().solve(global_system_load_vector);
	//global_displ = (global_system_stiffness_matrix.transpose() * global_system_stiffness_matrix).ldlt().solve(global_system_stiffness_matrix.transpose() * global_system_load_vector);
	//global_displ = global_system_stiffness_matrix.colPivHouseholderQr().solve(global_system_load_vector);
	//std::cout << "global displacement vector\n  " << global_displ << std::endl;

	/**
	* postprocessing
	* organize the generalized modal amplitude vector
	*/
	Eigen::MatrixXd V, V_x;
	amplitudeVx obVVx;
	obVVx.amplitude(mode_list, global_displ, LM_matrix, no_beam_elem, l, &V, &V_x);
	//std::cout << "V,x(x) matrix\n " << V_x << std::endl;
	//std::cout << "V(x) matrix\n " << V << std::endl;

	/**
	* calculation of displacements by combining the amplitude coefficients
	* coordinate transfer
	*/
	Eigen::MatrixXd u_x, v_y, w_z, u_X, v_Y, w_Z, X_disp, Y_disp, Z_disp;
	displacements obDis;
	obDis.amplitude(mode_list, V, V_x, no_beam_elem + 1, cross_sect_ref, r, &u_x, &v_y, &w_z, &u_X, &v_Y, &w_Z);
	//std::cout << "u_x\n " << u_x << std::endl;
	//std::cout << "v_y\n " << v_y << std::endl;
	//std::cout << "w_z\n " << w_z << std::endl;
	//std::cout << "u_X\n " << u_X << std::endl;
	//std::cout << "v_Y\n " << v_Y << std::endl;
	//std::cout << "w_Z\n " << w_Z << std::endl;
	X_disp = plot_scale * v_Y + x;
	Y_disp = plot_scale * u_X + y;
	Z_disp = plot_scale * w_Z + z;

	/**
	* Plot
	* displacement
	*/

	int nrow = X_disp.rows();
	int ncol = X_disp.cols();

	//double arr2[nrow][ncol];
	std::vector<std::vector<double>> X_data(nrow, std::vector<double>(ncol, 0)), Y_data(nrow, std::vector<double>(ncol, 0)), Z_data(nrow, std::vector<double>(ncol, 0));

#pragma omp parallel for 
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			X_data[i][j] = X_disp(i, j);
			Y_data[i][j] = Y_disp(i, j);
			Z_data[i][j] = Z_disp(i, j);

		}

	}
	//std::cout << "X \n " << X_data[nrow - 1][ncol - 1] << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	// Calculating total time taken by the program
	double time_taken =
		std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	time_taken *= 1e-9;
	std::cout << "Time taken by program is : " << std::fixed
		<< time_taken << std::setprecision(9);
	std::cout << " sec" << std::endl;
	
	matplot::surf(Y_data, Z_data, X_data)->face_alpha(0.5).edge_color("none");
	matplot::show();
	/*
    std::cout << "Do you want to run again (yes or n):" << std::endl;
	std::cin.ignore();
	std::cin.getline(again,256);
	std::cout << again << std::endl;
*/  
	} while (again == "yes"); 
}


