#pragma once

class GbtElementMatrix {
private:
	//Eigen::VectorXd loc_k1,loc_k2;
	Eigen::Matrix4d shLTxTxshL, shHTxxTxxshH, shHTTshH, shHTxTxshH, shHTxxTshH, shHTTxxshH;
	Eigen::MatrixXd element_matrix;
	int leng;

public:
	Eigen::MatrixXd ikLinear(Eigen::MatrixXd ik_tensor, double l,  int element_dof) {
		/*
		* Determines the linear element stiffness matrix,
		*/
		element_matrix.setZero(element_dof, element_dof);
		/*
		* The following stiffness matrices are built based on the Hermite shape functions 
        * The matrices are imported from MAPLE
        * The detail formulation of these metrices can be found on section 2.2.6 
		*/
		shLTxTxshL << 0.37e2 / 0.10e2 / l, -0.189e3 / 0.40e2 / l, 0.27e2 / 0.20e2 / l, -0.13e2 / 0.40e2 / l,
			-0.189e3 / 0.40e2 / l, 0.54e2 / 0.5e1 / l, -0.297e3 / 0.40e2 / l, 0.27e2 / 0.20e2 / l,
			0.27e2 / 0.20e2 / l, -0.297e3 / 0.40e2 / l, 0.54e2 / 0.5e1 / l, -0.189e3 / 0.40e2 / l,
			-0.13e2 / 0.40e2 / l, 0.27e2 / 0.20e2 / l, -0.189e3 / 0.40e2 / l, 0.37e2 / 0.10e2 / l;

		shHTxxTxxshH <<
			12 / pow(l, 3), 6 / pow(l, 2), -12 / pow(l, 3), 6 / pow(l, 2),
			6 / pow(l, 2), 4 / l, -6 / pow(l, 2), 2 / l,
			-12 / pow(l, 3), -6 / pow(l, 2), 12 / pow(l, 3), -6 / pow(l, 2),
			6 / pow(l, 2), 2 / l, -6 / pow(l, 2), 4 / l;

		shHTxTxshH <<
			0.6e1 / 0.5e1 / l, 0.1e1 / 0.10e2, -0.6e1 / 0.5e1 / l, 0.1e1 / 0.10e2,
			0.1e1 / 0.10e2, 0.2e1 / 0.15e2 * l, -0.1e1 / 0.10e2, -l / 30,
			-0.6e1 / 0.5e1 / l, -0.1e1 / 0.10e2, 0.6e1 / 0.5e1 / l, -0.1e1 / 0.10e2,
			0.1e1 / 0.10e2, -l / 30, -0.1e1 / 0.10e2, 0.2e1 / 0.15e2 * l;

		shHTTshH <<
			0.13e2 / 0.35e2 * l, 0.11e2 / 0.210e3 * pow(l,2), 0.9e1 / 0.70e2 * l, -0.13e2 / 0.420e3 * pow(l,2),
			0.11e2 / 0.210e3 * pow(l,2), pow(l,3) / 105, 0.13e2 / 0.420e3 * pow(l,2), -pow(l,3) / 140,
			0.9e1 / 0.70e2 * l, 0.13e2 / 0.420e3 * pow(l,2), 0.13e2 / 0.35e2 * l, -0.11e2 / 0.210e3 * pow(l,2),
			-0.13e2 / 0.420e3 * pow(l,2), -pow(l,3) / 140, -0.11e2 / 0.210e3 * pow(l,2), pow(l,3) / 105;

		shHTxxTshH <<
			-0.6e1 / 0.5e1 / l, -0.1e1 / 0.10e2, 0.6e1 / 0.5e1 / l, -0.1e1 / 0.10e2,
			-0.11e2 / 0.10e2, -0.2e1 / 0.15e2 * l, 0.1e1 / 0.10e2, l / 30,
			0.6e1 / 0.5e1 / l, 0.1e1 / 0.10e2, -0.6e1 / 0.5e1 / l, 0.1e1 / 0.10e2,
			-0.1e1 / 0.10e2, l / 30, 0.11e2 / 0.10e2, -0.2e1 / 0.15e2 * l;

		shHTTxxshH <<
			-0.6e1 / 0.5e1 / l, -0.11e2 / 0.10e2, 0.6e1 / 0.5e1 / l, -0.1e1 / 0.10e2,
			-0.1e1 / 0.10e2, -0.2e1 / 0.15e2 * l, 0.1e1 / 0.10e2, l / 30,
			0.6e1 / 0.5e1 / l, 0.1e1 / 0.10e2, -0.6e1 / 0.5e1 / l, 0.11e2 / 0.10e2,
			-0.1e1 / 0.10e2, l / 30, 0.1e1 / 0.10e2, -0.2e1 / 0.15e2 * l;

        leng = ik_tensor.rows();
#pragma omp parallel for 		
		for (int i = 0; i < leng ; i++) {
		//	loc_k1 = ik_tensor(i, Eigen::seq(6, 9));
		//	loc_k2 = ik_tensor(i, Eigen::seq(10, 13));
			if (ik_tensor(i, 0) == 1 && ik_tensor(i, 1) == 1) {
				element_matrix(ik_tensor(i, Eigen::seq(7, 10)), ik_tensor(i, Eigen::seq(11, 14))) += ik_tensor(i, 2) * shLTxTxshL;
			}
			if (ik_tensor(i, 0) != 1 && ik_tensor(i, 1) != 1) {
				element_matrix(ik_tensor(i, Eigen::seq(7, 10)), ik_tensor(i, Eigen::seq(11, 14))) += ik_tensor(i, 2) * shHTxxTxxshH + ik_tensor(i, 3) * shHTTshH + ik_tensor(i, 4) * shHTxTxshH + ik_tensor(i, 6) * shHTxxTshH + ik_tensor(i, 5) * shHTTxxshH;
			}
		};
		return element_matrix;
	};
};
