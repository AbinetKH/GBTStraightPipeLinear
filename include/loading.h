#pragma once

class GbtLoading {
private:
	int m, k_int, length_mode_list, loc_s, loc_e;
	std::string k, k_type, k_mode;
	Eigen::VectorXd element_load_vector;
	Eigen::Vector4d mode_load_vector;

public:
	Eigen::VectorXd loadVector(std::vector<std::string> mode_list, double r, double l, double q, Eigen::VectorXd LM_element, double ver_s, double ver_e) {
		/**
		* modal Decompostion
		*/
		length_mode_list = std::size(mode_list);
		element_load_vector.setZero(LM_element.rows());
		mode_load_vector << 6.0 * l, pow(l, 2), 6.0 * l, -pow(l, 2);
		mode_load_vector *= r*q / 12.0;

//#pragma omp parallel for 
		for (int i = 0; i < length_mode_list; i++) {
			loc_s = i * 4;
			loc_e = (i * 4) + 3;
			k = mode_list[i];
			k_type = k.substr(k.size() - 1, k.size());
			k_mode = k.substr(0, k.size() - 2);

			if (k_mode != "t" && k_mode != "a" && k_mode != "1") {
				k_int = std::stoi(k_mode);
				if (k_int % 2 == 0) {
					m = k_int / 2;
				}
				else if (k_int % 2 != 0) {
					m = (k_int - 1) / 2;
				}
			}
			if (k_type == "t" || k_type == "1" || k_type == "u" ) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) *= 0;
			}
			else if (k_type == "a") {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = ((-ver_s / 2.0 - 1.0 / 4.0 * sin(2.0 * ver_s)) - (-ver_e / 2.0 - 1.0 / 4.0 * sin(2.0 * ver_e))) * mode_load_vector;
			} 
			else if (k_int % 2 == 0 && k_type == "c" && m == 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = (-2.0 * pow(cos(ver_s), 4) - (-2.0 * pow(cos(ver_e), 4))) * mode_load_vector;
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) += 0.5 * element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0));
			}
			else if (k_int % 2 != 0 && k_type == "c" && m == 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = (( - 4.0 * (ver_s / 4.0 + 1.0 / 4.0 * sin(2.0 * ver_s) + 1.0 / 16.0 * sin(4.0 * ver_s))) -
					(-4.0 * (ver_e / 4.0 + 1.0 / 4.0 * sin(2.0 * ver_e) + 1.0 / 16.0 * sin(4.0 * ver_e))))* mode_load_vector;
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) += 0.5 * element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0));
			}
			else if (k_int % 2 == 0 && k_type == "c" && m != 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = ((pow(m, 2) * (-cos((-2.0 + m) * ver_s) / (4.0 * (-2.0 + m)) - cos(m * ver_s) / (2.0 * m) - cos((2.0 + m) * ver_s) / (4.0 * (2.0 + m)))) -
					(pow(m, 2) * (-cos((-2.0 + m) * ver_e) / (4.0 * (-2.0 + m)) - cos(m * ver_e) / (2.0 * m) - cos((2.0 + m) * ver_e) / (4.0 * (2.0 + m))))) * mode_load_vector;//z
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) += 0.5 * element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0));
			}
			else if (k_int % 2 != 0 && k_type == "c" && m != 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = ((-1.0 / 4.0 * pow(m, 2) * (sin((-2.0 + m) * ver_s) / (-2.0 + m) + (2.0 * sin(m * ver_s)) / m + sin((2.0 + m) * ver_s) / (2.0 + m))) -
					(-1.0 / 4.0 * pow(m, 2) * (sin((-2.0 + m) * ver_e) / (-2.0 + m) + (2.0 * sin(m * ver_e)) / m + sin((2.0 + m) * ver_e) / (2.0 + m)))) * mode_load_vector;
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) += 0.5 * element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0));
			}
			else if (k_int % 2 == 0 && k_type == "v" && m == 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = 0.5 * (-2.0 * pow(cos(ver_s), 4) - (-2.0 * pow(cos(ver_e), 4))) * mode_load_vector;
			}
			else if (k_int % 2 != 0 && k_type == "v" && m == 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = 0.5 * ((-4.0 * (ver_s / 4.0 + 1.0 / 4.0 * sin(2.0 * ver_s) + 1.0 / 16.0 * sin(4.0 * ver_s))) -
					(-4.0 * (ver_e / 4.0 + 1.0 / 4.0 * sin(2.0 * ver_e) + 1.0 / 16.0 * sin(4.0 * ver_e)))) * mode_load_vector;
			}
			else if (k_int % 2 == 0 && k_type == "v" && m != 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = 0.5 * ((pow(m, 2) * (-cos((-2.0 + m) * ver_s) / (4.0 * (-2.0 + m)) - cos(m * ver_s) / (2.0 * m) - cos((2.0 + m) * ver_s) / (4.0 * (2.0 + m)))) -
					(pow(m, 2) * (-cos((-2.0 + m) * ver_e) / (4.0 * (-2.0 + m)) - cos(m * ver_e) / (2.0 * m) - cos((2.0 + m) * ver_e) / (4.0 * (2.0 + m))))) * mode_load_vector;

			}
			else if (k_int % 2 != 0 && k_type == "v" && m != 2) {
				element_load_vector(LM_element(Eigen::seq(loc_s, loc_e), 0)) = 0.5 * ((-1.0 / 4.0 * pow(m, 2) * (sin((-2.0 + m) * ver_s) / (-2.0 + m) + (2.0 * sin(m * ver_s)) / m + sin((2.0 + m) * ver_s) / (2.0 + m))) -
					(-1.0 / 4.0 * pow(m, 2) * (sin((-2.0 + m) * ver_e) / (-2.0 + m) + (2.0 * sin(m * ver_e)) / m + sin((2.0 + m) * ver_e) / (2.0 + m)))) * mode_load_vector;
			}

		}

		return element_load_vector;
	}
	};