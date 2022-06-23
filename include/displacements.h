#pragma once
#include "GBT_function.h"
class displacements {
private:
    Eigen::MatrixXd loc_u_x, loc_v_y, loc_w_z, loc_u_X, loc_v_Y, loc_w_Z;
    Eigen::Matrix3d coordTrans, coordTrans_num;
    Eigen::Vector3d local_vec, global_vec;
    std::string k, k_mode;
    int length_mode_list, LM_len, LM_pos, LM_pos1, LM_pos2, LM_pos3;
    double vertheta, sum_u, sum_v, sum_w, gbt_func_u, gbt_func_v, gbt_func_w, varphi, pi = 3.14159, varphi_add=0;

public:
    Eigen::MatrixXd SphericalTransformationMatrix(double vertheta, double varphi) {
        coordTrans << -sin(vertheta) * sin(varphi), cos(varphi), cos(vertheta)* sin(varphi),
            -sin(vertheta)* cos(varphi), -sin(varphi), cos(vertheta)* cos(varphi),
            cos(vertheta), 0, sin(vertheta);
        return coordTrans;
    
    };
    void  amplitude(std::vector<std::string> mode_list, Eigen::MatrixXd V, Eigen::MatrixXd V_x, int no_of_nodes, int cross_sect_ref, double r, Eigen::MatrixXd* u_x, Eigen::MatrixXd* v_y, Eigen::MatrixXd* w_z, Eigen::MatrixXd* u_X, Eigen::MatrixXd* v_Y, Eigen::MatrixXd* w_Z) {
        length_mode_list = std::size(mode_list) ;
        varphi = 0.1 / (no_of_nodes-1);
        loc_u_x.setZero(no_of_nodes, cross_sect_ref);
        loc_v_y.setZero(no_of_nodes, cross_sect_ref);
        loc_w_z.setZero(no_of_nodes, cross_sect_ref);
        loc_u_X.setZero(no_of_nodes, cross_sect_ref);
        loc_v_Y.setZero(no_of_nodes, cross_sect_ref);
        loc_w_Z.setZero(no_of_nodes, cross_sect_ref);
        GbtFunction obGbtFunc;
        for (int i = 0; i < no_of_nodes; i++) {
            vertheta = 0;
            for (int j = 0; j < cross_sect_ref; j++) {
                sum_u = 0, sum_v = 0, sum_w = 0;
                for (int m = 0; m < length_mode_list; m++) {
                    k = mode_list[m];
                   // k_mode = k.substr(0, k.size() - 2);
                    gbt_func_u = obGbtFunc.deformationFuncU(k, r, vertheta);
                    gbt_func_v = obGbtFunc.deformationFuncV(k, r, vertheta);
                    gbt_func_w = obGbtFunc.deformationFuncW(k, r, vertheta);
                    sum_u += V_x(i, m) * gbt_func_u;
                    sum_v += V(i, m) * gbt_func_v;
                    sum_w += V(i, m) * gbt_func_w;
                }
                loc_u_x(i, j) = sum_u;
                loc_v_y(i, j) = sum_v;
                loc_w_z(i, j) = sum_w;
                local_vec << loc_v_y(i, j), loc_u_x(i, j), loc_w_z(i, j);
                coordTrans_num=SphericalTransformationMatrix(vertheta, varphi_add);
                global_vec = coordTrans_num * local_vec;
                loc_u_X(i, j) = global_vec(1);
                loc_v_Y(i, j) = global_vec(0);
                loc_w_Z(i, j) = global_vec(2);
                vertheta += (2 * pi / (cross_sect_ref - 1));
            }
            varphi_add += varphi;
        }
        *u_x = loc_u_x, * v_y = loc_v_y, * w_z = loc_w_z, * u_X = loc_u_X, * v_Y = loc_v_Y, * w_Z = loc_w_Z;
    };
};