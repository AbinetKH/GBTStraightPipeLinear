#pragma once

#include "mode_coupling_ik.h"

class couplingTensors {
private:
	int count = 0, LM_k1, LM_k2;
	double Cik, Bik, Dik, muDik_1, muDik_2;
	std::string k1, k2;
	Eigen::MatrixXd ik_tensor;
	GbtLinearModeCoupling obCoupling;
public:
	Eigen::MatrixXd ikTensor(std::vector<std::string> mode_list, Eigen::VectorXd LM_element, int length_mode_list, double mu, double t, double r, double Q, double K, double G) {
		/*
		* determines the second - order coupling tensor based on GbtLinearModeCoupling class
		*/ 
		ik_tensor.setZero(1, 15);
//#pragma omp parallel for 
		for (int i = 0; i < length_mode_list; i++) {
			for (int j = 0; j < length_mode_list; j++) {
				//std::cout << "modes all   " << mode_list[j] << std::endl;
				k1 = mode_list[i];
				k2 = mode_list[j];
				if (k1[0] == k2[0]) {
				   // std::cout << "coupled modes   " << k1.size()  << k2.substr(k2.size() - 1,k2.size()) << std::endl;
					Cik = obCoupling.ikCStraightPipe(k1, k2, t, r, Q, K, G);
					Bik = obCoupling.ikBStraightPipe(k1, k2, t, r, Q, K, G);
					Dik = obCoupling.ikDStraightPipe(k1, k2, t, r, Q, K, G);
					muDik_1 = mu * obCoupling.ikDmu1StraightPipe(k1, k2, t, r, Q, K, G);
					muDik_2 = mu * obCoupling.ikDmu2StraightPipe(k1, k2, t, r, Q, K, G);
					if (k1 == "1|1") {
						ik_tensor(count, 0) = 1;
					}
					if (k2 == "1|1") {
						ik_tensor(count, 1) = 1;
					}
					ik_tensor(count, 2) = Cik;
					ik_tensor(count, 3) = Bik;
					ik_tensor(count, 4) = Dik;
					ik_tensor(count, 5) = muDik_1;
					ik_tensor(count, 6) = muDik_2;
					LM_k1 = i * 4;
					LM_k2 = j * 4;
					// std::cout << "LM_element  " << LM_element(Eigen::seq(LM_k1, LM_k1 + 3), 0) << std::endl;
					ik_tensor(count, Eigen::seq(7, 10)) = LM_element(Eigen::seq(LM_k1, LM_k1 + 3), 0);
					ik_tensor(count, Eigen::seq(11, 14)) = LM_element(Eigen::seq(LM_k2, LM_k2 + 3), 0);
					ik_tensor.conservativeResize(ik_tensor.rows() + 1, Eigen::NoChange);
					ik_tensor(count + 1, Eigen::all) = 0 * ik_tensor(count, Eigen::all);
					count++; //task to do: make count independent for parallization
				}

			}
		}
		return ik_tensor;
	}
};


