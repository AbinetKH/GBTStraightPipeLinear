#pragma once

class GbtFunction {
private:
	int m, k_int;
	std::string k_type, k_mode;
public:
	double deformationFuncU(std::string k, double r, double vertheta) {
		/**
		* returns the GBT u(theta) function 
		*/
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

		if (k_type == "t" || k_type == "a" || k_type == "v") {
			return 0;
		}
		else if (k_type == "1" ) {
			return 1.0;
		}
		else if (k_int % 2 == 0 && (k_type == "c" || k_type == "u" )) {
			return (r * sin(m * vertheta));
		}
		else if (k_int % 2 != 0 && (k_type == "c" || k_type == "u")) {
			return (-r * cos(m * vertheta));
		}
	}
	double deformationFuncU1xdiff(std::string k, double r, double vertheta) {
		/**
        * returns the GBT u(theta) function first derivative
        */
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

		if (k_type == "t" || k_type == "a" || k_type == "v" || k_type == "1") {
			return 0;
		}
		else if (k_int % 2 == 0 && (k_type == "c" || k_type == "u")) {
			return (r * m * cos(m * vertheta));
		}
		else if (k_int % 2 != 0 && (k_type == "c" || k_type == "u")) {
			return (-r * m * sin(m * vertheta));
		}
	}
	double deformationFuncV(std::string k, double r, double vertheta) {
		/**
		* returns the GBT v(theta) function
		*/
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

		if (k_type == "t" ) {
			return r;
		}
		else if (k_type == "1" || k_type == "a" || k_type == "u") {
			return 0;
		}
		else if (k_int % 2 == 0 && (k_type == "c" || k_type == "v")) {
			return (-m * cos(m * vertheta));
		}
		else if (k_int % 2 != 0 && (k_type == "c" || k_type == "v")) {
			return (-m * sin(m * vertheta));
		}
	}
	double deformationFuncV1xdiff(std::string k, double r, double vertheta) {
		/**
		* returns the GBT v(theta) function first derivative
		*/
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

		if (k_type == "t" || k_type == "a" || k_type == "u" || k_type == "1") {
			return 0;
		}
		else if (k_int % 2 == 0 && (k_type == "c" || k_type == "v")) {
			return (pow(m, 2) * sin(m * vertheta));
		}
		else if (k_int % 2 != 0 && (k_type == "c" || k_type == "v")) {
			return (-pow(m, 2) * cos(m * vertheta));
		}
	}
	double deformationFuncW(std::string k, double r, double vertheta) {
		/**
		* returns the GBT w(theta) function
		*/
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


		if (k_type == "1" || k_type == "t" || k_type == "u" || k_type == "v") {
			return 0;
		}
		else if (k_type == "a") {
			return 1.0;
		}
		else if (k_int % 2 == 0 && k_type == "c" ) {
			return (-pow(m, 2) * sin(m * vertheta));
		}
		else if (k_int % 2 != 0 && k_type == "c" ) {
			return (pow(m, 2) * cos(m * vertheta)); 
		}
	}
	double deformationFuncW1xdiff(std::string k, double r, double vertheta) {
		/**
		* returns the GBT w(theta) function first derivative
		*/
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

		if (k_type == "1" || k_type == "t" || k_type == "u" || k_type == "v" || k_type == "a") {
			return 0;
		}
		else if (k_int % 2 == 0 && k_type == "c") {
			return (-pow(m, 3) * cos(m * vertheta));
		}
		else if (k_int % 2 != 0 && k_type == "c") {
			return (-pow(m, 3) * sin(m * vertheta));
		}
	}
	double deformationFuncW2xdiff(std::string k, double r, double vertheta) {
		/**
		* returns the GBT w(theta) function second derivative
		*/
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

		if (k_type == "1" || k_type == "t" || k_type == "u" || k_type == "v" || k_type == "a") {
			return 0;
		}
		else if (k_int % 2 == 0 && k_type == "c") {
			return (pow(m, 4) * sin(m * vertheta));
		}
		else if (k_int % 2 != 0 && k_type == "c") {
			return (-pow(m, 4) * cos(m * vertheta));
		}
	}
	double deformationFuncW3xdiff(std::string k, double r, double vertheta) {
		/**
		* returns the GBT w(theta) function third derivative
		*/
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

		if (k_type == "1" || k_type == "t" || k_type == "u" || k_type == "v" || k_type == "a") {
			return 0;
		}
		else if (k_int % 2 == 0 && k_type == "c") {
			return (pow(m, 5) * cos(m * vertheta));
		}
		else if (k_int % 2 != 0 && k_type == "c") {
			return (pow(m, 5) * sin(m * vertheta));
		}
	}
};