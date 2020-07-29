#include <iostream>

using namespace std;

#include "PQ.h"
#include "st.h"

int main() {
    PQ pq;
    /****************�򿪳�ʼ����ԭ�ļ�*****************/
    printf("Please Enter The Filename of The System:");
    gets(FILENAME);
    if ((fp = fopen(FILENAME, "r")) == NULL) {
        printf("Cannot Find The File:%s\n", FILENAME);
        printf("Press ENTER to Escape!");
        exit(0);
    } else {
        printf("success find file");
    }

    /********************��ϵͳ����������Ϣ*********************/
    fscanf(fp, "%d,%d,%d,%d,%lf,%lf", &N, &Nb, &Ng, &Nl, &V0, &epsilon);

    /*************��֧·������������ɡ�PV�ڵ�����ڴ�**********/
    Branch.resize((Nb + 1) * sizeof(struct Branch_Type));
    Generator.resize((Ng + 1) * sizeof(struct Generator_Type));
    Load.resize((Nl + 1) * sizeof(struct Load_Type));
    PVNode.resize(N * sizeof(struct PVNode_Type));

    pq.Datain(Nb, Nl, Ng, fp, Branch, Load, Generator, PVNode, Npv);  /*�ӳ�ʼ�ļ��ж�������*/

    for (n = 0; 1; n++)  /*��������������ļ�*/
    {
        if (FILENAME[n] == '.') {
            FILENAME[n] = '\0';
            strcat(FILENAME, "out.dat");
            break;
        }
    }
    if ((fp = fopen(FILENAME, "w")) == NULL)  /*�򿪼���������ļ�*/
    {
        printf("Cannot Find The File:%s\n", FILENAME);
        printf("Press ENTER to Escape!");
        exit(0);
    }

    /***********Ϊ���ɾ�������ڴ沢�γɵ��ɾ���***************/
    Yii.resize((N + 1) * sizeof(struct Yii_Type));
    Yiil.resize((N + 1) * sizeof(struct Yii_Type));
    Yij.resize((N + 1) * sizeof(struct Yij_Type));
    Yijl.resize((N + 1) * sizeof(struct Yij_Type));
    NYseq.resize((N + 1) * sizeof(int));
    pq.AdmittanceMatrix(N, Nb, Yii, Yiil, Yij, Yijl, Branch, NYseq);  /*���ɾ����γɺ���*/

    /*****�γ����ӱ�,����BX����(B'�к��Գ����ݺͷǱ�׼��ȣ�B"�к��Ե���)*****/
    U1.resize((N - 1) * (N - 2) / 2 * sizeof(struct U_Type));
    U2.resize((N - 1) * (N - 2) / 2 * sizeof(struct U_Type));
    D1.resize(N * sizeof(double));
    D2.resize(N * sizeof(double));
    NUsum1.resize(N * sizeof(int));
    NUsum2.resize(N * sizeof(int), 0);

    pq.Factorial(1, N, Npv, PVNode, NUsum1, Yii, Yij, NYseq, D1, U1);  /*�γ����ӱ�B'*/
    pq.AdmittanceMatrixAdd(Nb, Yii, Yiil, Yij, Yijl, Branch);  /*���ɾ���׷�ӽӵ�֧·����*/
    pq.Factorial(2, N, Npv, PVNode, NUsum2, Yiil, Yijl, NYseq, D2, U2);  /*�γ����ӱ�B"*/

    /****************������������õ����ӱ���е������******************/
    DI1.resize(N * sizeof(double));
    DI2.resize(N * sizeof(double));
    NodePow.resize((N + 1) * sizeof(struct NodalPow));
    NodeVol.resize((N + 1) * sizeof(struct NodalVol));

    /********************���͵�ѹ��ֵ***********************/

    for (i = 1; i <= N; i++) {
        NodeVol[i].V = V0;
        NodeVol[i].theta = 0.0;
    }
    for (n = 1; n <= Npv; n++) {
        i = PVNode[n].i;
        NodeVol[i].V = PVNode[n].V;
    }

    GenPower.resize((N + 1) * sizeof(struct GeneratorPower));
    fprintf(fp, "\t\t\tϵͳ���������� \n(1)�������̼�¼��\n��������\t    \t�й�����\t\t�޹�����\n\t\t   ��Pmax\tP-Node\t ��Qmax\t\tQ-Node\n");

    for (k = 0; 1; k++) {
        fprintf(fp, "  %2d:\t", k);
        if (Kp == 1) {
            pq.NodePower(1, N, NodeVol, NodePow, Yii, Yij, NYseq);  /*�ڵ㹦�ʼ��㺯��*/
            pq.Iteration(1, Generator, Load, PVNode, NodeVol, NodePow, GenPower, N, DI1, MaxError, ErrNode);  /*�������㺯��*/
            fprintf(fp, "\t%10.7lf\t %d\t", MaxError, ErrNode);
            if (MaxError >= epsilon)
                pq.FormulaSolution(1, U1, D1, NUsum1, DI1, N, NodeVol, V0);  /*���Է�������⺯��*/
            else
                Kp = 0;
        } else
            fprintf(fp, "\t\t\t\t");
        if (Kq == 1) {
            pq.NodePower(2, N, NodeVol, NodePow, Yii, Yij, NYseq);  /*�ڵ㹦�ʼ��㺯��*/
            pq.Iteration(2, Generator, Load, PVNode, NodeVol, NodePow, GenPower, N, DI2, MaxError, ErrNode);  /*�������㺯��*/
            fprintf(fp, "%10.7lf\t %d\n", MaxError, ErrNode);
            if (MaxError >= epsilon)
                pq.FormulaSolution(2, U2, D2, NUsum2, DI2, N, NodeVol, V0);  /*���Է�������⺯��*/
            else
                Kq = 0;
        } else
            fprintf(fp, "\n");
        if (Kp == 0 && Kq == 0)
            break;
        if (k > 1000) {
            fprintf(fp, "\n������������1000�Σ�ϵͳ������!\n");
        }
        printf("...");
    }

    fprintf(fp, "\n�ܵ�������Ϊ: %d ��!\n\n", k + 1);
    fprintf(fp, "(2)����������(�ڵ��ѹ)��\n Node\t\tV\t\t��\t\tP\t\tQ\n");
    pq.NodeDataOutput(fp, NodeVol, Generator, N, GenPower, NodePow, Load, Nl);
    fprintf(fp, "\n(3)����������(֧·����):\n    Branch\t\tP\t\tQ\n");
    pq.BranchDataOutput(fp, Nb, Branch, NodeVol);
    fprintf(fp, "\n\n\n");

    fclose(fp);
    return 0;
}
