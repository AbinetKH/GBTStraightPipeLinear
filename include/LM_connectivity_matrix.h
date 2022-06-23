#pragma once
#include <iostream>
#include <Eigen/Dense>


class LmConnectivityMatrix {
private:
	int length_mode_list, sum_b=0,sum_c,sum_tax,sum_high,sum_a;
	Eigen::MatrixXd LM_matrix;
public:
	Eigen::MatrixXd openEndPipe(std::vector<std::string> mode_list, int no_beam_elem) {
		/*
		* determines the LM connectivity matrix for assembly of the global stiffness matrix
		* refer Finite Element Procedures K. J. Bathe
		* code needs improvement
		*/
		length_mode_list = std::size(mode_list);
		LM_matrix.setZero(length_mode_list*4, no_beam_elem);
		
		for (int i = 0; i < no_beam_elem; i++){
			sum_c = 0, sum_tax=0, sum_high=0;
			for (int k = 0; k < length_mode_list; k++) {
				sum_a = 0;
				for (int j = 0; j < 4; j++) {
					//std::cout << "mode_listlm " << mode_list[k] << std::endl;
					if (mode_list[k] == "1|1") {
						LM_matrix(j + sum_c, i) = i + sum_a + sum_tax + sum_b + sum_high * (2 + no_beam_elem * 2);
						sum_a++;
					}
					else {
						LM_matrix(j + sum_c, i) = sum_high * (2 + no_beam_elem * 2) + sum_a + sum_b + sum_tax;
                        sum_a++;
					}
				}
				if (mode_list[k] == "1|1") {
					sum_tax +=((no_beam_elem * 4) - (no_beam_elem - 1));
				 }
				else{
					sum_high++ ;
				}
			    sum_c += 4;
			}
			sum_b += 2;
		}

		return LM_matrix;
	};
	Eigen::MatrixXd closeEndPipe(std::vector<std::string> mode_list, int no_beam_elem) {
		// to be completed, for toroidal pipe
		return LM_matrix;
	};
};