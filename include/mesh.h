#pragma once
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

class RectangularMesh {
private:
	Eigen::VectorXd r_d, w, r_d_cos, r_d_sin;
	Eigen::MatrixXd rr, cc;
	double  pi = 3.14159;
public:
	void meshStraightPipe(int cross_sect_ref, int no_beam_elem, double r, double length, Eigen::MatrixXd* y, Eigen::MatrixXd* z, Eigen::MatrixXd* x){
		r_d.setLinSpaced(cross_sect_ref, 0, 2 * pi);
		w.setLinSpaced(no_beam_elem+1, 0, length);
		rr.setOnes(no_beam_elem+1, 1);
		cc.setOnes(cross_sect_ref, 1);
		r_d_cos = r_d.array().cos() * r;
		r_d_sin = r_d.array().sin() * r;
		*y = rr * r_d_cos.transpose();
		*z = rr * r_d_sin.transpose();
		*x = w * cc.transpose();
	};
};

