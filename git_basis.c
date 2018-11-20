/****************************************************************************************************
 * linked-cell_basis_mod_cut
 * Marc Maaß, Juni 2018
 * 
 * Basis-Algorithmus für die Demonstration des Linked-Cell-Verfahren zur Bestimmung kurzreichweitiger Potentiale
 * Angewandt wird das Störmer-Verlet-Verfahren in Geschwindigkeitsform
 * 
 * * Zur Berechnung der Kräfte wird das Lennard-Jones-Potential genutzt. Entsprechend wurde die Funktion force2(...) angepasst
 * 
 * *_mod: die Funktion force(...) zur Berechnung der Kraftauswirkung eines Partikels j auf Partikel i wurde
 * angepasst um das 3. Newtonsche Gesetz auszunutzen (-> es müssen statt N^2, nur (N^2-N)/2 Operationen
 * ausgeführt werden. Die neue Funktion lautet force2(...)
 * 
 * *_cut: es wird der Abschneideradius eingeführt. Es genügt nur die Partikel, für die Kräfteberechnung eines 
 * Partikels, einzubeziehen, die einen wertigen Beitrag zur Kraft bzw. des Potentials beitragen. Partikel außerhalb 
 * eines definierten Radius
 *  
 * 
 * *************************************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//  Raumdimension
#define DIM 2
//  Definierung der Parameter sigma und epsilon für die Berechnung der Kräfte mit dem Lennard-Jones-Potential
#define sigma 1
#define epsilon 5
//  Abschneideradius r_cut definieren. Typischerweise rCut = 2.5 * epsilon
#define r_cut (2.5 * epsilon)
//  Makro fuer SQR
#define sqr(x) ((x) * (x))
//  datatype real
typedef double real;

typedef struct
{

    //  Masse
    real m;
    //  Position
    real x[DIM];
    //  Geschwindigkeit
    real v[DIM];
    //  Kraft
    real F[DIM];
    //  Kraft_old
    real F_old[DIM];
    //  Unterscheidung der Körper
    char body;

} Particle;

//  Routinen fuer einen Velocity-Stoermer-Verlet Zeitschritt fuer
//  ein EINZELNES Partikel
void updateX(Particle *p, real delta_t, real l[DIM])
{
    real a = delta_t * .5 / p->m;
    for (int d = 0; d < DIM; d++)
    {

        p->x[d] += delta_t * (p->v[d] + a * p->F[d]);
        
        //setze Randbedingung
        if(p->x[d] < 0)
            p->x[d] = 0;
        if(p->x[d] > l[d])
            p->x[d] = l[d];
        p->F_old[d] = p->F[d];
    }
}

void updateV(Particle *p, real delta_t)
{

    real a = delta_t * .5 / p->m;
    for (int d = 0; d < DIM; d++)
    {

        p->v[d] += a * (p->F[d] + p->F_old[d]);
    }
}

//  Routinen fuer einen Velocity-Stoermer-Verlet Zeitschritt fuer
//  einen VEKTOR von Partikeln
void compX_basis(Particle *p, int N, real delta_t, real l[DIM])
{

    for (int i = 0; i < N; i++)
    {

        updateX(&p[i], delta_t, l);
    }
}

void compV_basis(Particle *p, int N, real delta_t)
{

    for (int i = 0; i < N; i++)
    {

        updateV(&p[i], delta_t);
    }
}

void force2(Particle *i, Particle *j)
{

    real r = 0;
    for (int d = 0; d < DIM; d++)
    {

        r += sqr(j->x[d] - i->x[d]);
    }

    //  Wenn Partikel i, j einen größeren Abstand als r_cut aufweisen, kann die Kräfteberechnung zwischen den beiden
    //  vernachlässigt werde
    if (r <= r_cut)
    {
        real s = sqr(sigma) / r;
        s = sqr(s) * s;
        real f = 24 * epsilon * s / r * (1 - 2 * s);
        for (int d = 0; d < DIM; d++)
        {
            i->F[d] += f * (j->x[d] - i->x[d]);
            j->F[d] -= f * (j->x[d] - i->x[d]);
        }
    }
}

void compF_basis(Particle *p, int N)
{

    for (int i = 0; i < N; i++)
    {
        for (int d = 0; d < DIM; d++)
        {
            p[i].F[d] = 0;
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            force2(&p[i], &p[j]);
        }
    }
}

//  Berechnen der kinetischen Energie und printout
void compoutStatistic_basis(Particle *p, int N, real t)
{

    real e = 0;
    for (int i = 0; i < N; i++)
    {

        real v = 0;
        for (int d = 0; d < DIM; d++)
        {

            v += sqr(p[i].v[d]);
        }
        e += .5 * p[i].m * v;
    }

    printf("Kinetische Energie e zum Zeitpunkt t = %f: %f\n\n", t, e);
}

//  Darstellen des Zeitpunktes inkl. der Werte der Positionen und
//  Geschwindigkeiten der Partikel
void outputResults_basis(Particle *p, int N, real t)
{

    printf("Zeitpunkt t = %f:\n", t);
    for (int i = 0; i < N; i++)
    {
        printf("Particle No.%d:\t", i + 1);
        printf("m^%d= %f\t", i, p[i].m);
        for (int d = 0; d < DIM; d++)
        {
            printf("x^%d[%d]=\t%.3f\t", i + 1, d, p[i].x[d]);
        }

        for (int d = 0; d < DIM; d++)
        {

            printf("v^%d[%d]=\t%.3f\t", i + 1, d, p[i].v[d]);
        }
        printf("\n");
    }
}

void inputParameters_basis(real *delta_t, real *t_end, int *N_A, int *N_B, real *l1, real *l2)
{

    real delta_t_temp, t_end_temp;
    int N_A_temp, N_B_temp;
    real l_temp[DIM];

    printf("Insert delta_t:\n");
    scanf("%lf", &delta_t_temp);

    printf("Insert t_end:\n");
    scanf("%lf", &t_end_temp);

    printf("Insert Particlecount of body A: \n");
    scanf("%d", &N_A_temp);

    printf("Insert Particlecount of body B: \n");
    scanf("%d", &N_B_temp);

    printf("Insert Length of simulation domain in the x-dimension: \n");
    scanf("%lf", &l_temp[0]);

    printf("Insert Length of simulation domain in the y-dimension: \n");
    scanf("%lf", &l_temp[1]);

    printf("\n\t\t\t----------Parameter initialisiert----------\n\n");

    *delta_t = delta_t_temp;
    *t_end = t_end_temp;
    *N_A = N_A_temp;
    *N_B = N_B_temp;
    *l1 = l_temp[0];
    *l2 = l_temp[1];
}

//  Initialisierung der Partikeldaten
void initData_basis(Particle *p, int N_A, int N_B)
{

    int count = 0;

    //  Zugehörigkeit der Partikel zu Körper A oder B definieren
    for (int i = 0; i < N_A; i++)
    {

        p[i].body = 'A';
        p[i].m = 1;
        p[i].v[0] = 0.1;
        p[i].v[1] = 0.1;
    }

    for (int i = N_A; i < N_B + N_A; i++)
    {

        p[i].body = 'B';
        p[i].m = 1;
        p[i].v[0] = 0.1;
        p[i].v[1] = -15;
    }

    for (int i = 0; i < sqrt(N_A); i++)
    {

        for (int j = 0; j < sqrt(N_A); j++)
        {
            p[count].x[0] = i;
            p[count].x[1] = j;

            count++;
        }
    }

    for (int i = (sqrt(N_A) / 2) - (sqrt(N_B) / 2); i < (sqrt(N_A) / 2) - (sqrt(N_B) / 2) + sqrt(N_B); i++)
    {

        for (int j = sqrt(N_A) + 10; j < (sqrt(N_A) + 10) + sqrt(N_B); j++)
        {

            p[count].x[0] = i;
            p[count].x[1] = j;

            count++;
        }
    }
}

void outputCsv_basis(Particle *p, int N, int j)
{
    //Namen der CSV-Datei formatieren (Output.csv.TIMESTEP)
    char jChar[100], filename[100] = "outputtest/Output.csv.";
    itoa(j, jChar, 10);
    strcat(filename, jChar);
    printf("string lautet: %s", filename);

    //Filestream oeffnen
    FILE *fp;

    //Datei mit mit Bezeichnung filename erstellen/oeffnen
    fp = fopen(filename, "w+");

    fprintf(fp, "XCOORDS,YCOORDS,COLOR\n");

    for (int i = 0; i < N; i++)
    {
        if (p[i].body == 'A' || p[i].body == 'B')
        {

            fprintf(fp, "%f,%f,%f\n", p[i].x[0], p[i].x[1], round((p[i].v[0]+p[i].v[1])/2));
        }
        else
        {
            fprintf(fp, "%f,%f,%f\n", p[i].x[0], p[i].x[1],round(p[i].v[0]));
        }
    }
    fclose(fp);
}

void timeIntegration_basis(real t, real delta_t, real t_end, Particle *p, int N, char csv, real l[DIM])
{
    // Zaehlervariable fuer die Nummerierung des Outputs
    int i = 0;

    compF_basis(p, N);
    while (t < t_end)
    {

        t += delta_t;
        compX_basis(p, N, delta_t, l);
        compF_basis(p, N);
        compV_basis(p, N, delta_t);
        //compoutStatistic_basis(p, N, t);
        //outputResults_basis(p, N, t);
        if (csv == 'J' || csv == 'j')
        {
            outputCsv_basis(p, N, i);
        }
        else
        {
        }
        /*
        i++;
        if (t < t_end)
            printf("\n==================================================================================\n\n");
        else
            printf("\n\t\t\t----------Integration beendet----------\n\n");
        */
    }
}

void outputdata(real delta_t, real t_end, double timeused, int N_A, int N_B, real l[DIM]){

    FILE *fpdata;

    fpdata = fopen("results_basis.csv", "a");

    if(fpdata == NULL)
        printf("ERROR: Could not open file \"results_basis.csv\"");
    else{

        fprintf(fpdata, "\nBASIS,%f,%f,%d,%d,%.0f,%.0f,%.2f", delta_t,t_end,N_A,N_B,l[0],l[1],timeused);
        fclose(fpdata);
    }
}

int main()
{
    char csv;
    printf("Soll ein Output im .csv-Format generiert werden? [J/N]\n");
    scanf("%c", &csv);

    int N_A, N_B, N;
    real l[DIM];
    double timeused;

    real delta_t, t_end;
    inputParameters_basis(&delta_t, &t_end, &N_A, &N_B, &l[0], &l[1]);

    N = N_A + N_B;

    Particle *p = (Particle *)malloc(N * sizeof(*p));

    initData_basis(p, N_A, N_B);

    //Zeitmessung Beginn
    clock_t start = clock();
    timeIntegration_basis(0, delta_t, t_end, p, N, csv, l);
    clock_t end = clock();
    //Zeitmessung Ende

    timeused = (double)((end - start) / CLOCKS_PER_SEC);

    outputdata(delta_t, t_end, timeused, N_A, N_B, l);

    free(p);



    printf("\t\t\t\t   Genutzte Parameter:\n\n\t\tdelta_t = %f\tt_end = %.5f\tN = %d\t\n\n", delta_t, t_end, N);
    printf("\t\t -> Zeitschritte insgesamt: %.0f\n\n", (t_end / delta_t));
    printf("\t\t\t=======================================");
    printf("\n\n\t\t\tBenoetigte Zeit in Sekunden: %fs\n\n", timeused);
    printf("\t\t\t=======================================");

    printf("\t\t\t\tSUCCESS");

    return 0;
}