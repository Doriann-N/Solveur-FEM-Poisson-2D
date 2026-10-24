#pragma once
#include <iostream>
#include <vector>

struct Triplet{
    int i;
    int j;
    double value;
    Triplet(const int& p_i, const int& p_j, const double& p_value){
        i = p_i;
        j = p_j;
        value = p_value;    
    }
    bool operator<(const Triplet& other) const{
        if (i != other.i){
            return i < other.i;
        }
        return j < other.j;
    }
    void print() const{
        std::cout << "(i, j, value) : ( " << i << ", " << j << ", " << value <<" )\n";
    }
};

void display_mat(const std::vector<std::vector<double>> &K);

void display_vec(const std::vector<double> &v);

void display_vector_int(const std::vector<int>& v);

int ask_for_discretisation();
std::string ask_for_name_file();
bool ask_for_boolean(std::string what_to_ask);