#pragma once

class GbtLinearModeCoupling {
private:
	int m, k1_int, k2_int;
	double M_PI = 3.14159265358979323846;
	std::string k1_type, k2_type, k1_mode, k2_mode;
public:
	double ikCStraightPipe(std::string k1, std::string k2, double t, double r, double Q, double K, double G) {
		/**
		* determines the linear stiffness second-order coefficient tensor, ik_C
		* longitudinal stiffnesses coefficient with extensional (terms related with Q) and bending (terms related with K)
		* formulation of the equations used can be found on in Table 2.2 and 2.3
		*/
		k1_type = k1.substr(k1.size() - 1, k1.size());
		k2_type = k2.substr(k2.size() - 1, k2.size());
		k1_mode = k1.substr(0, k1.size() - 2);
		k2_mode = k2.substr(0, k2.size() - 2);
		if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1") {
			k1_int = std::stoi(k1_mode);
			if (k1_int % 2 == 0) {
				m = k1_int / 2;
			}
			else if (k1_int % 2 != 0) {
				m = (k1_int - 1) / 2;
			}
		}

		if (k1_mode == "a" && k1_mode == k2_mode) {
			return K*2*M_PI*r;
		}
		else if (k1_mode == "1" && k1_mode == k2_mode) {
			return Q * 2.0 * M_PI * r;
		}
		else if (k1_mode != "1" && k1_type == "c" && k2_type == "c" && k1_mode == k2_mode) {
			return (Q * M_PI * pow(r,3)) +(K*M_PI*r*pow(m,4));
		}
		else if (k1_mode != "1" && k1_type == "u" && k2_type == "u" && k1_mode == k2_mode) {
			return (Q * M_PI * pow(r, 3)) ;
		}
		else if (k1_mode != "1" && k1_type != "v" && k2_type != "v" && k1_mode == k2_mode && k1_type != k2_type) {
			return Q * M_PI * pow(r, 3);
		}
		else {
			return 0;
		}
	};
	double ikBStraightPipe(std::string k1, std::string k2, double t, double r, double Q, double K, double G) {
		/**
		* determines the linear stiffness second-order coefficient tensor, ik_B
		* transverse stiffnesses coefficient with extensional (terms related with Q) and bending (terms related with K)
		* formulation of the equations used can be found on in Table 2.2 and 2.3
		*/
		k1_type = k1.substr(k1.size() - 1, k1.size());
		k2_type = k2.substr(k2.size() - 1, k2.size());
		k1_mode = k1.substr(0, k1.size() - 2);
		k2_mode = k2.substr(0, k2.size() - 2);
		if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1") {
			k1_int = std::stoi(k1_mode);
			if (k1_int % 2 == 0) {
				m = k1_int / 2;
			}
			else if (k1_int % 2 != 0) {
				m = (k1_int - 1) / 2;
			}
		}

		if (k1_mode == "a" && k1_mode == k2_mode) {
			return Q * 2.0 * M_PI / r;
		}
		else if (k1_mode != "1" && k1_type == "c" && k2_type == "c" && k1_mode == k2_mode) {
			return (K * M_PI * pow(m, 4)/pow(r,3))*pow((pow(m,2)-1),2);
		}
		else if (k1_mode != "1" && k1_type == "v" && k2_type == "v" && k1_mode == k2_mode) {
			return  M_PI * pow(m, 4)*( Q/r+K/ pow(r, 3));
		}
		else if (k1_mode != "1" && k1_type != "u" && k2_type != "u" && k1_mode == k2_mode && k1_type != k2_type) {
			return (K * M_PI * pow(m, 4) / pow(r, 3)) * (1-pow(m, 2));
		}
		else {
			return 0;
		}
	};
	double ikDStraightPipe(std::string k1, std::string k2, double t, double r, double Q, double K, double G) {
		/**
		* determines the linear stiffness second-order coefficient tensor, ik_D
		* shear stiffnesses coefficient with extensional (terms related with Q) and bending (terms related with K)
		* formulation of the equations used can be found on in Table 2.2 and 2.3
		*/
		k1_type = k1.substr(k1.size() - 1, k1.size());
		k2_type = k2.substr(k2.size() - 1, k2.size());
		k1_mode = k1.substr(0, k1.size() - 2);
		k2_mode = k2.substr(0, k2.size() - 2);
		if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1") {
			k1_int = std::stoi(k1_mode);
			if (k1_int % 2 == 0) {
				m = k1_int / 2;
			}
			else if (k1_int % 2 != 0) {
				m = (k1_int - 1) / 2;
			}
		}

		if (k1_mode == "t" && k1_mode == k2_mode) {
			return G * M_PI * t * r * (2 * pow(r, 2) + (3.0 / 8.0) * pow(t, 2));
		}
		else if (k1_mode != "1" && k1_type == "c" && k2_type == "c" && k1_mode == k2_mode) {
			return (G * M_PI * pow(t, 3) * pow(m, 2) / (3.0 * r)) * pow((pow(m, 2) - 1), 2);
		}
		else if (k1_mode != "1" && k1_type == "v" && k2_type == "v" && k1_mode == k2_mode) {
			return  G * M_PI * t * pow(m, 2) * (r + (3.0 * pow(t, 2) / (16.0 * r)));
		}
		else if (k1_mode != "1" && k1_type == "u" && k2_type == "u" && k1_mode == k2_mode) {
			return  G * M_PI * t * pow(m, 2) * (r + (pow(t, 2) / (48.0 * r)));
		}
		else if (k1_mode != "1" && k1_type != "c" && k2_type != "c" && k1_mode == k2_mode && k1_type != k2_type) {
			return G * M_PI * t * pow(m, 2) * ((pow(t, 2) / (16.0 * r)) - r); 
		}
		else if (k1_mode != "1" && k1_type != "u" && k2_type != "u" && k1_mode == k2_mode && k1_type != k2_type) {
			return (G * M_PI * pow(t, 3) * pow(m, 2) / (4.0 * r)) * (1-pow(m, 2));
		}
		else if (k1_mode != "1" && k1_type != "v" && k2_type != "v" && k1_mode == k2_mode && k1_type != k2_type) {
			return (G * M_PI * pow(t, 3) * pow(m, 2) / (12.0 * r)) * (1 - pow(m, 2));
		}
		else {
			return 0;
		}
	};
	double ikDmu1StraightPipe(std::string k1, std::string k2, double t, double r, double Q, double K, double G) {
		/**
		* determines the linear stiffness second-order coefficient tensor, ik_D_mu
		* shear stiffnesses coefficient with extensional (terms related with Q) and bending (terms related with K)
		* formulation of the equations used can be found on in Table 2.2 and 2.3
		*/
		k1_type = k1.substr(k1.size() - 1, k1.size());
		k2_type = k2.substr(k2.size() - 1, k2.size());
		k1_mode = k1.substr(0, k1.size() - 2);
		k2_mode = k2.substr(0, k2.size() - 2);
		if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1") {
			k1_int = std::stoi(k1_mode);
			if (k1_int % 2 == 0) {
				m = k1_int / 2;
			}
			else if (k1_int % 2 != 0) {
				m = (k1_int - 1) / 2;
			}
		}

        if (k1_mode != "1" && k1_type == "c" && k2_type == "c" && k1_mode == k2_mode) {
			return (K * M_PI * pow(m, 4) / r) * (1 - pow(m, 2));
		}
		else if (k1_mode != "1" && k1_type == "v" && k2_type == "u" && k1_mode == k2_mode ) {
			return Q * M_PI * r * pow(m, 2) ;
		}
		else if (k1_mode != "1" && k1_type == "v" && k2_type == "c" && k1_mode == k2_mode ) {
			return M_PI * pow(m, 2) * (Q * r + (K * pow(m, 2)) / r);
		}
		else {
			return 0;
		}
	};
	double ikDmu2StraightPipe(std::string k1, std::string k2, double t, double r, double Q, double K, double G) {
		/**
		* determines the linear stiffness second-order coefficient tensor, ik_D_mu
		* shear stiffnesses coefficient with extensional (terms related with Q) and bending (terms related with K)
		* formulation of the equations used can be found on in Table 2.2 and 2.3
		*/
		k1_type = k1.substr(k1.size() - 1, k1.size());
		k2_type = k2.substr(k2.size() - 1, k2.size());
		k1_mode = k1.substr(0, k1.size() - 2);
		k2_mode = k2.substr(0, k2.size() - 2);
		if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1") {
			k1_int = std::stoi(k1_mode);
			if (k1_int % 2 == 0) {
				m = k1_int / 2;
			}
			else if (k1_int % 2 != 0) {
				m = (k1_int - 1) / 2;
			}
		}

		if (k1_mode != "1" && k1_type == "c" && k2_type == "c" && k1_mode == k2_mode) {
			return (K * M_PI * pow(m, 4) / r) * (1 - pow(m, 2));
		}
		else if (k1_mode != "1" && k1_type == "u" && k2_type == "v" && k1_mode == k2_mode) {
			return Q * M_PI * r * pow(m, 2);
		}
		else if (k1_mode != "1" && k1_type == "c" && k2_type == "v" && k1_mode == k2_mode) {
			return M_PI * pow(m, 2) * (Q * r + (K * pow(m, 2)) / r);
		}
		else {
			return 0;
		}
	};
};
