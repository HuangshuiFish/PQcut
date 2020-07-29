//
// Created by 17820 on 2020/7/29.
//

#ifndef PQ_CPP_PQ_H
#define PQ_CPP_PQ_H

#include <iostream>
#include <vector>
//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
#include <cstring>
#include <math.h>

using namespace std;

struct Branch_Type {
    int i, j;
    double R, X, YK;
};

struct Generator_Type {
    int i;
    double P, Q;
    double V;
};

struct Load_Type {
    int i;
    double P, Q;
    double V;
};

struct PVNode_Type {
    int i;
    double V;
};

struct Yii_Type {
    double G, B;
};

struct Yij_Type {
    double G, B;
    int j;
};

struct NodalPow {
    double P, Q;
};

struct NodalVol {
    double V, theta;
};

struct GeneratorPower {
    double P, Q;
};

struct U_Type {
    double value;
    int j;
};


class PQ {
public:
    PQ() = default;  // 明确告诉编译器，使用默认实现
    ~PQ() = default;  // 明确告诉编译器，使用默认实现
    void Datain(int Nb, int Nl, int Ng, FILE *fp, vector<Branch_Type> &Branch, vector<Load_Type> &Load,
                vector<Generator_Type> &Generator, vector<PVNode_Type> &PVNode, int &Npv);  /*从初始数据文件中读入数据*/
    void AdmittanceMatrix(int N, int Nb, vector<Yii_Type> &Yii, vector<Yii_Type> &Yiil, vector<Yij_Type> &Yij,
                          vector<Yij_Type> &Yijl, vector<Branch_Type> &Branch, vector<int> &NYseq1);  /*形成导纳矩阵函数*/
    void AdmittanceMatrixAdd(int Nb, vector<Yii_Type> &Yii, vector<Yii_Type> &Yiil, vector<Yij_Type> &Yij,
                             vector<Yij_Type> &Yijl, vector<Branch_Type> &Branch);  /*形成导纳矩阵追加接地支路函数*/
    void Factorial(int flag, int N, int Npv, vector<PVNode_Type> &PVNode, vector<int> &NUsum, vector<Yii_Type> &Yii,
                   vector<Yij_Type> &Yij, vector<int> &NYseq, vector<double> &D, vector<U_Type> &U);  /*形成因子表函数*/
    void NodePower(int flag, int N, vector<NodalVol> &NodeVol, vector<NodalPow> &NodePow, vector<Yii_Type> &Yii,
                   vector<Yij_Type> &Yij, vector<int> &NYseq);  /*节点功率计算函数*/
    void Iteration(int flag, vector<Generator_Type> &Generator, vector<Load_Type> &Load, vector<PVNode_Type> &PVNode,
                   vector<NodalVol> &NodeVol, vector<NodalPow> &NodePow, vector<GeneratorPower> &GenPower, int N,
                   vector<double> &DI, double &MaxError, int &ErrNode);  /*迭代计算函数*/
    void FormulaSolution(int flag, vector<U_Type> &U, vector<double> &D, vector<int> &NUsum, vector<double> &DI, int N,
                         vector<NodalVol> &NodeVol, double V0);  /*线性方程组求解函数*/
    void NodeDataOutput(FILE *fp, vector<NodalVol> &NodeVol, vector<Generator_Type> &Generator, int N,
                        vector<GeneratorPower> &GenPower, vector<NodalPow> &NodePow, vector<Load_Type> &Load,
                        int Nl);  /*节点数据输出函数*/
    void BranchDataOutput(FILE *fp, int Nb, vector<Branch_Type> &Branch, vector<NodalVol> &NodeVol);  /*支路数据输出函数*/
};


#endif //PQ_CPP_PQ_H
