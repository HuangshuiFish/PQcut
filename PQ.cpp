//
// Created by 17820 on 2020/7/29.
//

#include "PQ.h"


void PQ::Datain(int Nb, int Nl, int Ng, FILE *fp, vector<Branch_Type> &Branch, vector<Load_Type> &Load,
                vector<Generator_Type> &Generator, vector<PVNode_Type> &PVNode, int &Npv) {
    int n;
    /*********************从初始数据文件中读取支路、发电机、负荷数据************************/
    for (n = 1; n <= Nb; n++) {
        fscanf(fp, "%d,%d,%lf,%lf,%lf", &(Branch[n].i), &(Branch[n].j), &(Branch[n].R), &(Branch[n].X),
               &(Branch[n].YK));
    }
    for (n = 1; n <= Ng; n++) {
        fscanf(fp, "%d,%lf,%lf,%lf", &(Generator[n].i), &(Generator[n].P), &(Generator[n].Q), &(Generator[n].V));
        if ((Generator[n].V) < 0) {
            (Npv)++;
            PVNode[(Npv)].i = Generator[n].i;
            PVNode[(Npv)].V = -(Generator[n].V);
        }
    }
    for (n = 1; n <= Nl; n++) {
        fscanf(fp, "%d,%lf,%lf,%lf", &(Load[n].i), &(Load[n].P), &(Load[n].Q), &(Load[n].V));
        if ((Load[n].V) < 0) {
            (Npv)++;
            PVNode[(Npv)].i = Load[n].i;
            PVNode[(Npv)].V = -(Load[n].V);
        }
    }
    fclose(fp);
}


void PQ::AdmittanceMatrix(int N, int Nb, vector<Yii_Type> &Yii, vector<Yii_Type> &Yiil, vector<Yij_Type> &Yij,
                          vector<Yij_Type> &Yijl, vector<Branch_Type> &Branch, vector<int> &NYseq1) {
    /*****************不考虑接地支路时，形成节点导纳矩阵******************/
    int i, n;
//    int *NYsum = static_cast<int *>(malloc((N + 1) * sizeof(int)));
//    int *NYseq = static_cast<int *>(malloc((N + 1) * sizeof(int)));
    vector<int> NYsum;
    NYsum.resize((N + 1) * sizeof(int));
    vector<int> NYseq;
    NYseq.resize((N + 1) * sizeof(int));
    for (i = 1; i <= N; i++) {
        Yii[i].G = 0.0;
        Yii[i].B = 0.0;
        Yiil[i].G = 0.0;
        Yiil[i].B = 0.0;
        NYsum[i] = 0;
    }
    for (n = 1; n <= Nb; n++) {
        int i = abs(Branch[n].i);
        int j = abs(Branch[n].j);
        double R = (Branch[n].R);
        double X = (Branch[n].X);
        double YK = (Branch[n].YK);

        double Zmag2 = R * R + X * X;
        double Gij = R / Zmag2;
        double Bij = -X / Zmag2;

        double b_ij = -1.0 / X;

        if ((Branch[n].i < 0) || (Branch[n].j < 0)) {
            Yij[n].G = -Gij / YK;
            Yij[n].B = -Bij / YK;
            Yijl[n].G = 0.0;
            Yijl[n].B = -b_ij / YK;
        } else {
            Yij[n].G = -Gij;
            Yij[n].B = -Bij;
            Yijl[n].G = 0.0;
            Yijl[n].B = -b_ij;
        }
        Yij[n].j = j;
        Yijl[n].j = j;

        if ((Branch[n].i < 0) || (Branch[n].j < 0)) {
            Yii[i].G = Yii[i].G + Gij / YK;
            Yii[i].B = Yii[i].B + Bij / YK;
            Yii[j].G = Yii[j].G + Gij / YK;
            Yii[j].B = Yii[j].B + Bij / YK;

            Yiil[i].B = Yiil[i].B + b_ij / YK;
            Yiil[j].B = Yiil[j].B + b_ij / YK;
        } else {
            Yii[i].G = Yii[i].G + Gij;
            Yii[i].B = Yii[i].B + Bij;
            Yii[j].G = Yii[j].G + Gij;
            Yii[j].B = Yii[j].B + Bij;

            Yiil[i].B = Yiil[i].B + b_ij;
            Yiil[j].B = Yiil[j].B + b_ij;
        }
        NYsum[i] += 1;
    }

    NYseq[1] = 1;
    for (i = 1; i <= N - 1; i++) {
        NYseq[i + 1] = NYseq[i] + NYsum[i];
    }
    for (n = 1; n <= N; n++) {
        NYseq1[n] = NYseq[n];
    }

}

void PQ::AdmittanceMatrixAdd(int Nb, vector<Yii_Type> &Yii, vector<Yii_Type> &Yiil, vector<Yij_Type> &Yij,
                             vector<Yij_Type> &Yijl,
                             vector<Branch_Type> &Branch) {
    int n;
    /***********************在导纳矩阵中追加接地支路************************/
    for (n = 1; n <= Nb; n++) {
        int i = Branch[n].i;
        int j = Branch[n].j;
        double YK = Branch[n].YK;

        if (!(i < 0 || j < 0)) {
            double Bij = YK / 2.0;
            double b_ij = YK / 2.0;

            Yii[i].B = Yii[i].B + Bij;
            Yii[j].B = Yii[j].B + Bij;

            Yiil[i].B = Yiil[i].B + b_ij;
            Yiil[j].B = Yiil[j].B + b_ij;
        } else {
            double Gij;
            double Bij;
            double b_ij;
            if (i < 0) {
                i = abs(i);
                Gij = Yij[n].G;
                Bij = Yij[n].B;

                b_ij = Yijl[n].B;

                Yii[i].G = Yii[i].G + (1.0 - 1.0 / YK) * Gij;
                Yii[i].B = Yii[i].B + (1.0 - 1.0 / YK) * Bij;
                Yiil[i].B = Yiil[i].B + (1.0 - 1.0 / YK) * b_ij;


                Yii[j].G = Yii[j].G + (1.0 - YK) * Gij;
                Yii[j].B = Yii[j].B + (1.0 - YK) * Bij;
                Yiil[j].B = Yiil[j].B + (1.0 - YK) * b_ij;
            } else {
                j = abs(j);
                Gij = Yij[n].G;
                Bij = Yij[n].B;

                b_ij = Yijl[n].B;

                Yii[j].G = Yii[j].G + (1.0 - 1.0 / YK) * Gij;
                Yii[j].B = Yii[j].B + (1.0 - 1.0 / YK) * Bij;
                Yiil[j].B = Yiil[j].B + (1.0 - 1.0 / YK) * b_ij;


                Yii[i].G = Yii[i].G + (1.0 - YK) * Gij;
                Yii[i].B = Yii[i].B + (1.0 - YK) * Bij;
                Yiil[i].B = Yiil[i].B + (1.0 - YK) * b_ij;
            }
        }
    }

}

void PQ::Factorial(int flag, int N, int Npv, vector<PVNode_Type> &PVNode, vector<int> &NUsum, vector<Yii_Type> &Yii,
                   vector<Yij_Type> &Yij,
                   vector<int> &NYseq, vector<double> &D, vector<U_Type> &U) {

    /**************************因子表形成函数***************************/
    int n_pv, i_pv, j, n_u, i_above;
    int i, count;
    vector<double> B;
    B.resize((N + 1) * sizeof(double));
    double Btemp;
    n_pv = 1;
    i_pv = PVNode[1].i;

    for (i = 1; i < N; i++) {
        if ((flag == 2) && (i == i_pv)) {
            n_pv++;
            i_pv = PVNode[n_pv].i;
            NUsum[i] = 0;
            D[i] = 0.0;
            continue;
        } else {
            for (count = i + 1; count < N; count++)
                B[count] = 0.0;

            B[i] = Yii[i].B;
            for (count = NYseq[i]; count < NYseq[i + 1]; count++) {
                j = Yij[count].j;
                B[j] = Yij[count].B;
            }
            if (flag == 2) {
                for (count = 1; count <= Npv; count++) {
                    j = PVNode[count].i;
                    B[j] = 0.0;
                }

            }


            n_u = 1;
            i_above = 1;
            while (i_above <= (i - 1)) {
                count = 1;

                while (count <= NUsum[i_above]) {
                    if (U[n_u].j == i) {

                        Btemp = U[n_u].value / D[i_above];
                        while (count <= NUsum[i_above]) {
                            j = U[n_u].j;
                            B[j] = B[j] - Btemp * U[n_u].value;
                            count++;
                            n_u++;
                        }
                        break;
                    }
                    count++;
                    n_u++;
                }
                i_above++;
            }
            Btemp = 1.0 / B[i];
            D[i] = Btemp;

            count = 0;
            for (j = i + 1; j < N; j++) {
                if (B[j] != 0.0) {
                    U[n_u].value = B[j] * Btemp;
                    U[n_u].j = j;
                    count++;
                    n_u++;
                }
            }
            NUsum[i] = count;
        }
    }
}

void PQ::NodePower(int flag, int N, vector<NodalVol> &NodeVol, vector<NodalPow> &NodePow, vector<Yii_Type> &Yii,
                   vector<Yij_Type> &Yij, vector<int> &NYseq) {
/*********************计算各节点功率函数**********************/
    double A, B, Vi;
    int i, n, j;
    double VV, theta;
    for (i = 1; i <= N; i++) {
        if (flag == 1)
            NodePow[i].P = 0.0;
        else
            NodePow[i].Q = 0.0;
    }

    for (i = 1; i <= N; i++) {
        Vi = NodeVol[i].V;

        if (flag == 1) {
            A = Yii[i].G;
        } else {
            A = -Yii[i].B;
        }

        if (flag == 1)
            NodePow[i].P += Vi * Vi * A;
        else { NodePow[i].Q += Vi * Vi * A; }

        if (i == N) {
            break;
        } else {
            for (n = NYseq[i]; n <= NYseq[i + 1] - 1; n++) {
                if (flag == 1) {
                    A = Yij[n].G;
                    B = Yij[n].B;
                } else {
                    A = -Yij[n].B;
                    B = Yij[n].G;
                }
                j = Yij[n].j;

                VV = Vi * NodeVol[j].V;
                theta = NodeVol[i].theta - NodeVol[j].theta;
                A = A * VV * cos(theta);
                B = B * VV * sin(theta);
                if (flag == 1) {
                    NodePow[i].P += (A + B);
                    NodePow[j].P += (A - B);
                } else {
                    NodePow[i].Q += (A + B);
                    NodePow[j].Q += (A - B);
                }
            }
        }
    }
}

void PQ::Iteration(int flag, vector<Generator_Type> &Generator, vector<Load_Type> &Load, vector<PVNode_Type> &PVNode,
                   vector<NodalVol> &NodeVol, vector<NodalPow> &NodePow, vector<GeneratorPower> &GenPower, int N,
                   vector<double> &DI,
                   double &MaxError, int &ErrNode) {
    /*****************迭代计算函数*****************/
    int i = 1, n_g = 1, n_l = 1, n_pv = 1, i_g = Generator[1].i, i_l = Load[1].i, i_pv = PVNode[1].i;
    double Vi, Wi, Wtemp;
    (MaxError) = 0.0;

    do {
        Vi = NodeVol[i].V;
        if (i == i_l) {
            if (flag == 1) {
                Wi = Load[n_l].P;
            } else {
                Wi = Load[n_l].Q;
            }
            n_l += 1;
            i_l = Load[n_l].i;
        } else {
            Wi = 0.0;
        }
        Wtemp = Wi;
        if (flag == 1)
            Wi = Wi - NodePow[i].P;
        else
            Wi = Wi - NodePow[i].Q;

        if (i == i_g) {
            if (flag == 1) {
                NodePow[i].P = Wtemp;
                GenPower[i_g].P = -Wi;
            } else {
                NodePow[i].Q = Wtemp;
                GenPower[i_g].Q = -Wi;
            }

            if (flag == 1) {
                Wi += Generator[n_g].P;
            } else {
                Wi += Generator[n_g].Q;
            }
            n_g += 1;
            i_g = Generator[n_g].i;
        }

        if (i == N) {
            break;
        } else {
            if (flag == 2 && i == i_pv) {
                n_pv += 1;
                i_pv = PVNode[n_pv].i;
                DI[i] = 0.0;
            } else {
                if (fabs(Wi) > (MaxError)) {
                    (MaxError) = fabs(Wi);
                    (ErrNode) = i;
                }

                DI[i] = Wi / Vi;
            }
        }
        i += 1;
    } while (1);
}

void PQ::FormulaSolution(int flag, vector<U_Type> &U, vector<double> &D, vector<int> &NUsum, vector<double> &DI, int N,
                         vector<NodalVol> &NodeVol,
                         double V0) {
    /******************进行线性方程组的求解********************/
    int n_u;
    int i, count;
    int j;
    double DItemp, Dtheta, DV;
    n_u = 1;

    for (i = 1; i <= N - 1; i++) {
        DItemp = DI[i];
        for (count = 1; count <= NUsum[i]; count++) {
            j = U[n_u].j;
            DI[j] = DI[j] - DItemp * U[n_u].value;
            n_u++;
        }
        DI[i] = DItemp * D[i];
    }

    for (i = N - 1; i >= 1; i--) {
        DItemp = DI[i];
        for (count = 1; count <= NUsum[i]; count++) {
            n_u -= 1;
            j = U[n_u].j;
            DItemp = DItemp - DI[j] * U[n_u].value;
        }
        DI[i] = DItemp;
    }

    for (i = 1; i <= N - 1; i++) {
        if (flag == 1) {
            Dtheta = DI[i] / V0;
            NodeVol[i].theta = NodeVol[i].theta - Dtheta;
        } else {
            DV = DI[i];
            NodeVol[i].V -= DV;
        }
    }
}

void PQ::NodeDataOutput(FILE *fp, vector<NodalVol> &NodeVol, vector<Generator_Type> &Generator, int N,
                        vector<GeneratorPower> &GenPower, vector<NodalPow> &NodePow, vector<Load_Type> &Load, int Nl) {
    /***************节点数据输出**********************/
    double Vmin = NodeVol[1].V;
    double V, theta, P, Q;
    int i_g = Generator[1].i;
    int VminNode = 1;
    int n_g = 1;
    int i;

    for (i = 1; i <= N; i++) {
        theta = NodeVol[i].theta / 3.14159 * 180;
        V = NodeVol[i].V;
        if (V < Vmin) {
            Vmin = V;
            VminNode = i;
        }
        if (i == i_g) {
            P = GenPower[i].P;
            Q = GenPower[i].Q;
            n_g += 1;
            i_g = Generator[n_g].i;
        } else {
            P = 0.0;
            Q = 0.0;
        }
        if (i != N)
            fprintf(fp, "  %d\t   %10.7lf\t   %10.7lf\t  %10.7lf\t  %10.7lf\n", i, V, theta, P, Q);
        else
            fprintf(fp, "  %d\t   %10.7lf\t   %10.7lf\t  %10.7lf\t  %10.7lf\n", i, V, theta, NodePow[i].P - Load[Nl].P,
                    NodePow[i].Q - Load[Nl].Q);

    }
    fprintf(fp, "系统最低电压=%10.7lf,节点=%d\n", Vmin, VminNode);
}

void PQ::BranchDataOutput(FILE *fp, int Nb, vector<Branch_Type> &Branch, vector<NodalVol> &NodeVol) {
    /***************支路数据输出***************/
    double PLoss = 0.0, QLoss = 0.0;
    int n;
    int i, j;
    double R, X, YK, Y, theta, Ei, Ej, Fi, Fj, Vi, Vj, DE, DF;
    double Zmag2, Ir, Ii;
    double Pij, Qij, Pji, Qji;
    for (n = 1; n <= Nb; n++) {
        i = abs(Branch[n].i);
        j = abs(Branch[n].j);
        R = Branch[n].R;
        X = Branch[n].X;
        YK = Branch[n].YK;

        Vi = NodeVol[i].V;
        theta = NodeVol[i].theta;
        Ei = Vi * cos(theta);
        Fi = Vi * sin(theta);

        Vj = NodeVol[j].V;
        theta = NodeVol[j].theta;
        Ej = Vj * cos(theta);
        Fj = Vj * sin(theta);

        if (Branch[n].i < 0 || Branch[n].j < 0) {
            if (Branch[n].i < 0) {
                Ei = Ei / YK;
                Fi = Fi / YK;
            } else {
                Ej = Ej / YK;
                Fj = Fj / YK;
            }
            YK = 0.0;
        }

        DE = Ei - Ej;
        DF = Fi - Fj;
        Zmag2 = R * R + X * X;
        Ir = (DE * R + DF * X) / Zmag2;
        Ii = (DF * R - DE * X) / Zmag2;

        Pij = Ir * Ei + Ii * Fi;
        Qij = Ir * Fi - Ii * Ei;

        Pji = -(Ir * Ej + Ii * Fj);
        Qji = -(Ir * Fj - Ii * Ej);

        Qij -= (Vi * Vi * YK / 2.0);
        Qji -= (Vj * Vj * YK / 2.0);

        PLoss = PLoss + Pij + Pji;
        QLoss = QLoss + Qij + Qji;

        fprintf(fp, "  %3d->%3d\t   %10.7lf\t   %10.7lf\n  %3d->%3d\t   %10.7lf\t   %10.7lf\n", i, j, Pij, Qij, j, i,
                Pji, Qji);
    }

    fprintf(fp, "     损耗\t   %10.7lf\t   %10.7lf\n", PLoss, QLoss);
}
