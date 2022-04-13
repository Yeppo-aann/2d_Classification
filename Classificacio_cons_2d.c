//Treball C Anna Danot Sanchez
//NIU: 1514857

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIM 3
#define TOL 1e-6
#define PI 3.1415926


double producteescalar(double vect_A[], double vect_B[]);
void productevectorial (double vect_A[], double vect_B[], double vect_C[]);
double norma(double vect[]);
double determinant(double vect_A[], double vect_B[], double vect_C[]);
double angle(double vect_A[], double vect_B[], double vect_C[]);
int perpendicular(double vect_A[], double vect_B[], double vect_C[]);
void imprimeixintv (int l, unsigned int vect[]);
int tipus2d(int n, double Mv[], double u[], unsigned int vmin[]);


int main ()
{
    double *Mv;
    int files,columnes;
    char nomf[20];
    int n;
    int c=0;
    unsigned int vmin[4];
    double u[DIM]; //vector perpendicular a tots

    FILE * fitxer;

    printf("Cal entrar el nom del fitxer acabat amb .dat\n");
    printf("Entra el nom del fitxer: ");
    scanf("%s",nomf);

    fitxer = fopen(nomf, "r");

    if (fitxer==NULL)
    {
        printf("No s'ha pogut obrir el fitxer %s",nomf);
        return -1;
    }
    int e;
    while (EOF!=(e=fgetc(fitxer)))
    {
        if (e ==',')
            c++;
    }
    n=c/2;  //Llegeix en nombre de vectors en funci¢ de les comes que hi han al fitxer
    printf("\nS'ha llegit un fitxer amb %d vectors\n",n);
    fclose(fitxer);

    fitxer = fopen(nomf, "r");
    files=DIM;
    columnes=n;
    Mv = (double*) malloc(files*columnes*sizeof(double));

    int f=0;
    while (EOF!=fscanf(fitxer, "%lf, %lf, %lf", &Mv[f], &Mv[n+f], &Mv[2*n+f]))
    {
        f++;
    }
    fclose(fitxer);

    int k=0;
    double k1,k2,k3;
    double v1[]={Mv[0],Mv[n],Mv[2*n]};
    for (int i=0;i<n;i++)
    {
        if (k==0)
        {
            double X[]={Mv[i],Mv[n+i],Mv[2*n+i]};
            k1=v1[0]/X[0];
            k2=v1[1]/X[1];
            k3=v1[2]/X[2];
            if (k1!=k2)
                k=i;
            if (k1!=k3)
                k=i;
            if (k2!=k3)
                k=i;
        }
    }
    double K[]={Mv[k],Mv[n+k],Mv[2*n+k]};
    productevectorial(v1,K,u);
    //Tot i que el vector perpendicular u dongui nul, aix¢ vol dir que tots els vectors son linealment dependents i per tant
    //ser… o el con nul, o una recta, o una semirecta, i en aquests casos no importa el signe de l'angle, per tant deixem u com vector nul.

    for (int i=0; i<n; i++)
    {
        double T[DIM]={Mv[i], Mv[n+i], Mv[2*n+i]};
        if (perpendicular(T,u,u)==0)
        {
            printf("No %cs de tipus 2d",130);
            return -1;
        }
    }

    int resultat;
    resultat= tipus2d(n, Mv, u, vmin);
    if (resultat == 0)
        printf("\nEl con generat %cs el {0}\n",130);
    if (resultat == 1)
        printf("\nEl con generat %cs una semirecta\n",130);
    if (resultat == 2)
        printf("\nEl con generat %cs una recta\n",130);
    if (resultat == 3)
        printf("\nEl con generat %cs un angle pla\n",130);
    if (resultat == 4)
        printf("\nEl con generat %cs un semipla\n",130);
    if (resultat == 5)
        printf("\nEl con generat %cs un pla\n",130);
    return 0;

}



double producteescalar(double vect_A[], double vect_B[])
{
    int i;
    double producte=0;
    for (i=0; i<DIM; i++)
        producte = producte + vect_A[i]*vect_B[i];
    return producte;
}

void productevectorial (double vect_A[], double vect_B[], double vect_C[DIM])
{
    vect_C[0]=vect_A[1]*vect_B[2]-vect_A[2]*vect_B[1];
    vect_C[1]=vect_A[2]*vect_B[0]-vect_A[0]*vect_B[2];
    vect_C[2]=vect_A[0]*vect_B[1]-vect_A[1]*vect_B[0];
    return;
}

double norma(double vect[])
{
    return sqrt(producteescalar(vect,vect));
}

double determinant(double vect_A[], double vect_B[], double vect_C[])
{
    int i;
    double M[DIM][DIM];
    for (i=0; i<DIM; i++)
    {
        M[0][i]=vect_A[i];
        M[1][i]=vect_B[i];
        M[2][i]=vect_C[i];
    }
    return (M[0][0]*M[1][1]*M[2][2]+M[2][0]*M[0][1]*M[1][2]+M[1][0]*M[0][2]*M[2][1]-(M[2][0]*M[1][1]*M[0][2]+M[2][1]*M[0][0]*M[1][2]+M[0][0]*M[2][1]*M[1][2]));
}

double angle(double vect_A[], double vect_B[], double vect_C[])
{
    double x,x1,a;
    x=producteescalar(vect_A,vect_B)/(norma(vect_A)*norma(vect_B));
    x1=x*1e13;
    a=acos(trunc(x1)/1e13);
    if (determinant(vect_A,vect_B,vect_C)!=0)
        return a*(determinant(vect_A,vect_B,vect_C)/fabs(determinant(vect_A,vect_B,vect_C)));
    else
        return fabs(a);
}

int perpendicular(double vect_A[], double vect_B[],double vect_C[])
{
    if (PI/2-TOL < fabs(angle(vect_A,vect_B,vect_C)))
        if (fabs(angle(vect_A,vect_B,vect_C)) < PI/2+TOL)
        return 1;
    if(fabs(angle(vect_A,vect_B,vect_C)) < PI/2-TOL || fabs(angle(vect_A,vect_B,vect_C)) > PI/2+TOL)
        return 0;
    return -1;
}

void imprimeixintv(int l, unsigned int vect[])
{
    int i;
    printf("Un sistema minimal de generadors %cs: ",130);
    for(i=0;i<l+1;i++)
    {
         printf(" v_%d ",vect[i]+1);
    }
    return;
}

int tipus2d( int n, double Mv[], double u[], unsigned int vmin[])
{
    int i, j;
    int dr, esq;
    double Mangle[n][n];
    double anglemax;
    int lvmin=0;

    if (n==0)
        return 0;
    if (n==1)
    {
        vmin[lvmin]=0;
        imprimeixintv(lvmin,vmin);
        FILE * generadors;
        generadors = fopen("sortida.dat","w");
        if (generadors == NULL)
        {
            printf("No s'ha pogut crear el fitxer..\n\n");
            return -1;
        }
        for (i=0;i<lvmin+1; i++)
            fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
        fclose(generadors);
        return 1; //Semirecta
    }
    if (n>1)
    {
        dr=0;
        esq=1;
        double V0[DIM]={Mv[0],Mv[n],Mv[2*n]};
        double V1[DIM]={Mv[1],Mv[n+1],Mv[2*n+1]};
        anglemax=angle(V0, V1, u);
        for (i=0; i<n; i++)
        {
            for (j=0; j<n; j++)
            {
                double U[DIM]={Mv[i], Mv[n+i], Mv[2*n+i]};
                double V[DIM]={Mv[j], Mv[n+j], Mv[2*n+j]};
                Mangle[i][j]=angle(U, V, u);
                if (Mangle[i][j]>anglemax)
                {
                    dr=i;
                    esq=j;
                    anglemax=Mangle[i][j];
                }
            }
        }

        if (anglemax>0-TOL && anglemax<0+TOL) //anglemax=0
        {
            vmin[lvmin]=0;
            imprimeixintv(lvmin,vmin);
            FILE * generadors;
            generadors = fopen("sortida.dat","w");
            if (generadors == NULL)
            {
                printf("No s'ha pogut crear el fitxer..\n\n");
                return -1;
            }
            for (i=0;i<lvmin+1; i++)
                fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
            fclose(generadors);
            return 1; //Semirecta
        }
        if (anglemax<PI) //pla o anglepla
        {
            for (i=0; i<n; i++)
            {
                if (Mangle[dr][i]<0)
                {
                    printf("He aribat aqui 1");
                    vmin[lvmin]=dr;
                    lvmin++;
                    vmin[lvmin]=esq;
                    lvmin++;
                    vmin[lvmin]=i;
                    imprimeixintv(lvmin,vmin);
                    FILE * generadors;
                    generadors = fopen("sortida.dat","w");
                    if (generadors == NULL)
                    {
                        printf("No s'ha pogut crear el fitxer..\n\n");
                        return -1;
                    }
                    for (i=0;i<lvmin+1; i++)
                        fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
                    fclose(generadors);
                    return 5;//pla
                }
            }
            vmin[lvmin]=dr;
            lvmin++;
            vmin[lvmin]=esq;
            imprimeixintv(lvmin,vmin);
            FILE * generadors;
            generadors = fopen("sortida.dat","w");
            if (generadors == NULL)
            {
                printf("No s'ha pogut crear el fitxer..\n\n");
                return -1;
            }
            for (i=0;i<lvmin+1; i++)
                fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
            fclose(generadors);
            return 3;//anglepla
        }
        else //=PI
        {
            for (i=0; i<n; i++)
            {
                for (j=0; j<n; j++)
                {
                    if (Mangle[i][j]<PI)
                    {
                        if (Mangle[i][j]<0-TOL || Mangle[i][j]>0+TOL  )
                        {
                            for (int k=0; k<n; k++)
                            {
                                if (Mangle[dr][k]/fabs(Mangle[dr][k])!= Mangle[i][j]/fabs(Mangle[i][j]))
                                {
                                    if(0-TOL>Mangle[dr][k] || Mangle[dr][k]>0+TOL)
                                    {
                                        vmin[lvmin]=dr;
                                        lvmin++;
                                        vmin[lvmin]=esq;
                                        lvmin++;
                                        vmin[lvmin]=j;
                                        lvmin++;
                                        vmin[lvmin]=k;
                                        imprimeixintv(lvmin,vmin);
                                        FILE * generadors;
                                        generadors = fopen("sortida.dat","w");
                                        if (generadors == NULL)
                                        {
                                            printf("No s'ha pogut crear el fitxer..\n\n");
                                            return -1;
                                        }
                                        for (i=0;i<lvmin+1; i++)
                                            fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
                                        fclose(generadors);
                                        return 5;//pla
                                    }
                                }
                            }
                            vmin[lvmin]=dr;
                            lvmin++;
                            vmin[lvmin]=esq;
                            lvmin++;
                            vmin[lvmin]=j;
                            imprimeixintv(lvmin,vmin);
                            FILE * generadors;
                            generadors = fopen("sortida.dat","w");
                            if (generadors == NULL)
                            {
                                printf("No s'ha pogut crear el fitxer..\n\n");
                                return -1;
                            }
                            for (i=0;i<lvmin+1; i++)
                                fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
                            fclose(generadors);
                            return 4;//Semipla
                        }
                    }

                }
            }
            vmin[lvmin]=dr;
            lvmin++;
            vmin[lvmin]=esq;
            imprimeixintv(lvmin,vmin);
            FILE * generadors;
            generadors = fopen("sortida.dat","w");
            if (generadors == NULL)
            {
                printf("No s'ha pogut crear el fitxer..\n\n");
                return -1;
            }
            for (i=0;i<lvmin+1; i++)
                fprintf(generadors, "v(%d)=(%lf,%lf,%lf)\n",i+1,Mv[vmin[i]],Mv[n+vmin[i]],Mv[2*n+vmin[i]]);
            fclose(generadors);
            return 2; //Recta
        }
    }

    if (n<0)
    {
        return -1;
    }
    return -1;
}


