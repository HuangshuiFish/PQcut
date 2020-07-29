#include <iostream>

using namespace std;

#include "PQ.h"
#include "st.h"

int main() {
    PQ pq;
    /****************打开初始数据原文件*****************/
    printf("Please Enter The Filename of The System:");
    gets(FILENAME);
    if ((fp = fopen(FILENAME, "r")) == NULL) {
        printf("Cannot Find The File:%s\n", FILENAME);
        printf("Press ENTER to Escape!");
        exit(0);
    } else {
        printf("success find file");
    }

    /********************读系统基本数据信息*********************/
    fscanf(fp, "%d,%d,%d,%d,%lf,%lf", &N, &Nb, &Ng, &Nl, &V0, &epsilon);

    /*************给支路、发电机、负荷、PV节点分配内存**********/
    Branch.resize((Nb + 1) * sizeof(struct Branch_Type));
    Generator.resize((Ng + 1) * sizeof(struct Generator_Type));
    Load.resize((Nl + 1) * sizeof(struct Load_Type));
    PVNode.resize(N * sizeof(struct PVNode_Type));

    pq.Datain(Nb, Nl, Ng, fp, Branch, Load, Generator, PVNode, Npv);  /*从初始文件中读入数据*/

    for (n = 0; 1; n++)  /*建立计算结果输出文件*/
    {
        if (FILENAME[n] == '.') {
            FILENAME[n] = '\0';
            strcat(FILENAME, "out.dat");
            break;
        }
    }
    if ((fp = fopen(FILENAME, "w")) == NULL)  /*打开计算结果输出文件*/
    {
        printf("Cannot Find The File:%s\n", FILENAME);
        printf("Press ENTER to Escape!");
        exit(0);
    }

    /***********为导纳矩阵分配内存并形成导纳矩阵***************/
    Yii.resize((N + 1) * sizeof(struct Yii_Type));
    Yiil.resize((N + 1) * sizeof(struct Yii_Type));
    Yij.resize((N + 1) * sizeof(struct Yij_Type));
    Yijl.resize((N + 1) * sizeof(struct Yij_Type));
    NYseq.resize((N + 1) * sizeof(int));
    pq.AdmittanceMatrix(N, Nb, Yii, Yiil, Yij, Yijl, Branch, NYseq);  /*导纳矩阵形成函数*/

    /*****形成因子表,采用BX方案(B'中忽略充电电容和非标准变比，B"中忽略电阻)*****/
    U1.resize((N - 1) * (N - 2) / 2 * sizeof(struct U_Type));
    U2.resize((N - 1) * (N - 2) / 2 * sizeof(struct U_Type));
    D1.resize(N * sizeof(double));
    D2.resize(N * sizeof(double));
    NUsum1.resize(N * sizeof(int));
    NUsum2.resize(N * sizeof(int), 0);

    pq.Factorial(1, N, Npv, PVNode, NUsum1, Yii, Yij, NYseq, D1, U1);  /*形成因子表B'*/
    pq.AdmittanceMatrixAdd(Nb, Yii, Yiil, Yij, Yijl, Branch);  /*导纳矩阵追加接地支路函数*/
    pq.Factorial(2, N, Npv, PVNode, NUsum2, Yiil, Yijl, NYseq, D2, U2);  /*形成因子表B"*/

    /****************下面利用所求得的因子表进行迭代求解******************/
    DI1.resize(N * sizeof(double));
    DI2.resize(N * sizeof(double));
    NodePow.resize((N + 1) * sizeof(struct NodalPow));
    NodeVol.resize((N + 1) * sizeof(struct NodalVol));

    /********************先送电压初值***********************/

    for (i = 1; i <= N; i++) {
        NodeVol[i].V = V0;
        NodeVol[i].theta = 0.0;
    }
    for (n = 1; n <= Npv; n++) {
        i = PVNode[n].i;
        NodeVol[i].V = PVNode[n].V;
    }

    GenPower.resize((N + 1) * sizeof(struct GeneratorPower));
    fprintf(fp, "\t\t\t系统潮流计算结果 \n(1)迭代过程纪录：\n迭代次数\t    \t有功迭代\t\t无功迭代\n\t\t   ΔPmax\tP-Node\t ΔQmax\t\tQ-Node\n");

    for (k = 0; 1; k++) {
        fprintf(fp, "  %2d:\t", k);
        if (Kp == 1) {
            pq.NodePower(1, N, NodeVol, NodePow, Yii, Yij, NYseq);  /*节点功率计算函数*/
            pq.Iteration(1, Generator, Load, PVNode, NodeVol, NodePow, GenPower, N, DI1, MaxError, ErrNode);  /*迭代计算函数*/
            fprintf(fp, "\t%10.7lf\t %d\t", MaxError, ErrNode);
            if (MaxError >= epsilon)
                pq.FormulaSolution(1, U1, D1, NUsum1, DI1, N, NodeVol, V0);  /*线性方程组求解函数*/
            else
                Kp = 0;
        } else
            fprintf(fp, "\t\t\t\t");
        if (Kq == 1) {
            pq.NodePower(2, N, NodeVol, NodePow, Yii, Yij, NYseq);  /*节点功率计算函数*/
            pq.Iteration(2, Generator, Load, PVNode, NodeVol, NodePow, GenPower, N, DI2, MaxError, ErrNode);  /*迭代计算函数*/
            fprintf(fp, "%10.7lf\t %d\n", MaxError, ErrNode);
            if (MaxError >= epsilon)
                pq.FormulaSolution(2, U2, D2, NUsum2, DI2, N, NodeVol, V0);  /*线性方程组求解函数*/
            else
                Kq = 0;
        } else
            fprintf(fp, "\n");
        if (Kp == 0 && Kq == 0)
            break;
        if (k > 1000) {
            fprintf(fp, "\n迭代次数超过1000次，系统不收敛!\n");
        }
        printf("...");
    }

    fprintf(fp, "\n总迭代次数为: %d 次!\n\n", k + 1);
    fprintf(fp, "(2)潮流计算结果(节点电压)：\n Node\t\tV\t\tθ\t\tP\t\tQ\n");
    pq.NodeDataOutput(fp, NodeVol, Generator, N, GenPower, NodePow, Load, Nl);
    fprintf(fp, "\n(3)潮流计算结果(支路功率):\n    Branch\t\tP\t\tQ\n");
    pq.BranchDataOutput(fp, Nb, Branch, NodeVol);
    fprintf(fp, "\n\n\n");

    fclose(fp);
    return 0;
}
