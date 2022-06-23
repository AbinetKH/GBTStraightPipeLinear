#pragma once

class BoundaryCondition {
private:
    Eigen::VectorXd fixed_dofs;
    std::string k1, k1_mode;
    int count=0,count_dof=0;

public:
    Eigen::VectorXd dofsFixed(std::vector<std::string> mode_list, int length_mode_list , int no_beam_elem, std::vector<int> torsion_dof, std::vector<int> axisymmetric_dof, std::vector<int> extension_dof, std::vector<int> bending_dof, std::vector<int> allocal_dof, Eigen::MatrixXd LM_matrix) {
        for (int i = 0; i < length_mode_list; i++) {
            k1 = mode_list[i];
            k1_mode = k1.substr(0, k1.size() - 2);
            if (k1 == "t|t") {
                if (torsion_dof[0] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof-1)) = LM_matrix(Eigen::seq(count, count + 1), 0);
                }
                if (torsion_dof[1] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof - 1)) = LM_matrix(Eigen::seq(count + 2, count + 3), no_beam_elem - 1);
                }
                count += 4;
            }
            else if (k1 == "a|a") {
                if (axisymmetric_dof[0] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof-1)) = LM_matrix(Eigen::seq(count, count + 1), 0);
                }
                if (axisymmetric_dof[1] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof-1)) = LM_matrix(Eigen::seq(count + 2, count + 3), no_beam_elem - 1);
                }
                count += 4;
            }
            else if (k1 == "1|1") {
                if (extension_dof[0] == 1) {
                    count_dof += 1;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(count_dof-1) = LM_matrix(count, 0);
                }
                if (extension_dof[1] == 1) {
                    count_dof += 1;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(count_dof-1) = LM_matrix(count + 3, no_beam_elem - 1);
                }
                count += 4;
            }
            else if (k1_mode == "2" || k1_mode == "3") {
                if (bending_dof[0] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof - 1)) = LM_matrix(Eigen::seq(count, count + 1), 0);
                }
                if (bending_dof[1] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof - 1)) = LM_matrix(Eigen::seq(count + 2, count + 3), no_beam_elem - 1);
                }
                count += 4;
            }
            else if (k1_mode != "t" && k1_mode != "a" && k1_mode != "1" && k1_mode != "2" && k1_mode != "3") {
                if (allocal_dof[0] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof - 1)) = LM_matrix(Eigen::seq(count, count + 1), 0);
                }
                if (allocal_dof[1] == 1) {
                    count_dof += 2;
                    fixed_dofs.conservativeResize(count_dof);
                    fixed_dofs(Eigen::seq(count_dof - 2, count_dof - 1)) = LM_matrix(Eigen::seq(count + 2, count + 3), no_beam_elem - 1);
                }
                count += 4;
            }

        }
        return fixed_dofs;
    }
};

