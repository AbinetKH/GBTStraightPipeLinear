#pragma once
class amplitudeVx {
private:
    Eigen::MatrixXd loc_V, loc_V_x;
    std::string k, k_mode;
    int length_mode_list, LM_len, LM_pos, LM_pos1, LM_pos2, LM_pos3;

public:
    void  amplitude(std::vector<std::string> mode_list, Eigen::VectorXd global_displ, Eigen::MatrixXd LM_matrix, int no_beam_elem, int l, Eigen::MatrixXd* V, Eigen::MatrixXd* V_x) {
        length_mode_list = std::size(mode_list);
        loc_V_x.setZero(no_beam_elem + 1, length_mode_list);
        loc_V.setZero(no_beam_elem + 1, length_mode_list);
        LM_len = LM_matrix.cols()-1;
//#pragma omp parallel for 
        for (int i = 0; i < length_mode_list; i++) {
            k = mode_list[i];
            LM_pos = i * 4;
            LM_pos1 = LM_pos + 1;
            LM_pos2 = LM_pos + 2;
            LM_pos3 = LM_pos + 3;
            k_mode = k.substr(0, k.size() - 2);
            if (k_mode == "1") {
                loc_V_x(Eigen::all, i) = global_displ(Eigen::seq(LM_matrix(LM_pos, 0), LM_matrix(LM_pos3, LM_len), 3));
                            }
            else {               
                loc_V_x(Eigen::all, i) = global_displ(Eigen::seq(LM_matrix(LM_pos1, 0), LM_matrix(LM_pos3, LM_len), 2));
                loc_V(Eigen::all, i) = global_displ(Eigen::seq(LM_matrix(LM_pos, 0), LM_matrix(LM_pos2, LM_len), 2));
            }
        }
        *V = loc_V;
        *V_x = loc_V_x;
        }
    };