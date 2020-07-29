//
// Created by 17820 on 2020/7/29.
//

#ifndef PQ_CPP_ST_H
#define PQ_CPP_ST_H
int N;                             /*系统节点总数*/
int Nb;                            /*系统中支路数*/
int Ng;                            /*发电机节点总数,不包括平衡机*/
int Nl;                            /*负荷节点总数*/
double V0;                         /*系统标称电压*/
double epsilon;                    /*迭代收敛所要求的精确度*/
vector<Branch_Type> Branch;/*支路数据数组*/
vector<Generator_Type> Generator;  /*发电机数据数组*/
vector<Load_Type> Load;           /*负荷数据数组*/
vector<PVNode_Type> PVNode;        /*PV节点数据数组*/
int Npv = 0;                         /*PV节点个数*/
vector<Yii_Type> Yii;        /*导纳对角矩阵数组*/
vector<Yii_Type> Yiil;
vector<Yij_Type> Yij;        /*导纳非对角矩阵数组*/
vector<Yij_Type> Yijl;
vector<int> NYseq;                        /*矩阵非对角元素记录数组*/
vector<NodalPow> NodePow;          /*节点功率数组*/
vector<NodalVol> NodeVol;          /*节点电压数组*/
vector<GeneratorPower> GenPower;  /*发电机功率数组*/
vector<double> DI1;  /*存放误差项的数组*/
vector<double> DI2;
double MaxError = 0.0;               /*每次迭代过程中最大误差项*/
int ErrNode;                       /*最大误差项所在节点*/
int Kp = 1, Kq = 1;                     /*有功迭代和无功迭代的标记*/
int n, k, count;                     /*K为迭代次数*/
int i;
vector<U_Type> U1;            /*存放因子表非对角元素的数组*/
vector<U_Type> U2;
vector<double> D1;                    /*存放因子表对角元素的数组*/
vector<double> D2;
vector<int> NUsum1;               /*矩阵非对角元素记录数组*/
vector<int> NUsum2;

char FILENAME[20];
FILE *fp;
#endif //PQ_CPP_ST_H
