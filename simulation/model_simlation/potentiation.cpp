#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>

using namespace std;

# define gr (double)1.0 

void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2);
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2);
double t_min_fun(double pf, double pr, double *para, double *para2);
double v_min(double t0, double t, double pf0, double pr0, double *para, double *para2);
int compute_ht(double T_on, double T_off, int m, double *p0, double *para, double *para2);
double compute_rt(double T_on, double T_off, int m, double *p0, double *para, double *para2);
int compute_potentiation(double T_on, double T_off, double *para, double *para2);

int main(int argc, const char* argv[])
{
    double T_on = 8., T_off;
    double T_off_lis[3] = {8., 16., 24.};
    int T_off_lis_size = 3;
    int ht;
    /*
    para: 
    para[0]=gf, para[1]=gfs, para[2]=ade,
    para[3]=ades, para[4]=Dif, para[5]=gam, para[6]=rv
    para2:
    para2[0]=Pf_u, para2[1]=Pr_u,
    para2[2]=Pf_s, para2[3]=Pr_s,
    para2[4]=lam1, para2[5]=lam2,
    para2[6]=v11, para2[7]=v12, v1=(v11, v12)^t
    para2[8]=v21, para2[9]=v22, v2=(v21, v22)^t
    */
    double para[7]={0.0010012365164467779, 11.973267355969925, 0.022442256672103103, 0.28663730935606846, 0.49092280681082606, 0.005781634243102351, 0.07087592120432085};
    double para2[10]={((para[2]+para[6]*para[4])*para[0]+para[4]*gr)/((para[2]+para[4]+para[6]*para[4])*para[2]), 
                    (para[6]*para[4]*para[0] + (para[2]+para[4])*gr)/((para[2]+para[4]+para[6]*para[4])*para[2]), 
                    ((para[2]+para[6]*para[4])*para[1]+para[4]*gr)/((para[2]+para[6]*para[4])*para[3]+para[2]*para[4]), 
                    (para[6]*para[4]*para[1] + (para[3]+para[4])*gr)/((para[2]+para[6]*para[4])*para[3]+para[2]*para[4]), 
                    -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])-0.5*sqrt((para[3]-para[2]+para[4]-para[6]*para[4])*(para[3]-para[2]+para[4]-para[6]*para[4])+4.*para[6]*para[4]*para[4]), 
                    -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])+0.5*sqrt((para[3]-para[2]+para[4]-para[6]*para[4])*(para[3]-para[2]+para[4]-para[6]*para[4])+4.*para[6]*para[4]*para[4]), 
                    -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])-0.5*sqrt((para[3]-para[2]+para[4]-para[6]*para[4])*(para[3]-para[2]+para[4]-para[6]*para[4])+4.*para[6]*para[4]*para[4])+para[2]+para[6]*para[4], 
                    para[6]*para[4], 
                    para[4], 
                    -0.5*(para[3]+para[2]+para[4]+para[6]*para[4])+0.5*sqrt((para[3]-para[2]+para[4]-para[6]*para[4])*(para[3]-para[2]+para[4]-para[6]*para[4])+4.*para[6]*para[4]*para[4])+para[3]+para[4]};

    int m = 2;
    double *p0 = (double*)malloc(m*sizeof(double *));
    if(p0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    p0[0] = para2[0];
    p0[1] = para2[1];
    printf("\n#### Potentiation of habituation ####\n\n");
    for (int i = 0; i < T_off_lis_size; i++)
    {   
        T_off = T_off_lis[i];
        printf("T_s=%f, T=%f\n", T_on, T_on+T_off);
        ht = compute_ht(T_on, T_off, m, p0, para, para2);
        printf("ht=%d (first habituation)\n", ht);
        ht = compute_potentiation(T_on, T_off, para, para2);
        printf("ht=%d (second habituation)\n\n", ht);
    }

    free(p0);

    return 0;
}

// analytical solution for h(t)=0
void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2)
{  
    x[0] = para2[0]+(x0[0]-x0[1]-para2[0]+para2[1])*exp(-(para[2]+para[4]+para[6]*para[4])*(t-t0))/(1.+para[6]) + (para[6]*x0[0]+x0[1]-para[6]*para2[0]-para2[1])*exp(-para[2]*(t-t0))/(1.+para[6]);
    x[1] = para2[1]-para[6]*(x0[0]-x0[1]-para2[0]+para2[1])*exp(-(para[2]+para[4]+para[6]*para[4])*(t-t0))/(1.+para[6]) + (para[6]*x0[0]+x0[1]-para[6]*para2[0]-para2[1])*exp(-para[2]*(t-t0))/(1.+para[6]);
}

// analytical solution for h(t)=1
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2)
{
    x[0] = para2[2] + para2[6]*(para2[9]*(x0[0]-para2[2])-para2[8]*(x0[1]-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[8]*(-para2[7]*(x0[0]-para2[2])+para2[6]*(x0[1]-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);
    x[1] = para2[3] + para2[7]*(para2[9]*(x0[0]-para2[2])-para2[8]*(x0[1]-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[9]*(-para2[7]*(x0[0]-para2[2])+para2[6]*(x0[1]-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);
}

double t_min_fun(double pf, double pr, double *para, double *para2)
{
    /*
    para: 
    para[0]=gf, para[1]=gfs, para[2]=ade,
    para[3]=ades, para[4]=Dif, para[5]=gam, para[6]=rv
    para2:
    para2[0]=Pf_u, para2[1]=Pr_u,
    para2[2]=Pf_s, para2[3]=Pr_s,
    para2[4]=lam1, para2[5]=lam2,
    para2[6]=v11, para2[7]=v12, v1=(v11, v12)^t
    para2[8]=v21, para2[9]=v22, v2=(v21, v22)^t
    */
   double D = (para[3]-para[2]+para[4]-para[6]*para[4])*(para[3]-para[2]+para[4]-para[6]*para[4])+4.*para[6]*para[4]*para[4];
   double x1 = para2[4]*(para2[9]*(pf-para2[2])-para2[8]*(pr-para2[3]))*(para2[4]+para[2]);
   double x2 = para2[5]*(-para2[7]*(pf-para2[2])+para2[6]*(pr-para2[3]))*(para2[5]+para[3]);

   return log(x1/x2)/sqrt(D);
}

// velocity min
double v_min(double t0, double t, double pf0, double pr0, double *para, double *para2)
{
    double pf, pr;
    pf = para2[2] + para2[6]*(para2[9]*(pf0-para2[2])-para2[8]*(pr0-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[8]*(-para2[7]*(pf0-para2[2])+para2[6]*(pr0-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);
    pr = para2[3] + para2[7]*(para2[9]*(pf0-para2[2])-para2[8]*(pr0-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[9]*(-para2[7]*(pf0-para2[2])+para2[6]*(pr0-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);

    return para[5]*(pr-pf);
}

int compute_ht(double T_on, double T_off, int m, double *p0, double *para, double *para2)
{
    double t_min, res_1, res_2, v_base; 
    int ht;

    double *x0 = (double*)malloc(m*sizeof(double));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }
    double *x = (double*)malloc(m*sizeof(double));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    ht = 0;
    for (int i = 0; i < m; i++)
    {
        x0[i] = p0[i];
    }
    v_base = para[5]*(para2[1]-para2[0]);
    t_min = t_min_fun(x0[0], x0[1], para, para2);
    res_1 = v_base-v_min(0., t_min, x0[0], x0[1], para, para2);
    sol_1(0., T_on, x0, x, para, para2);
    sol_0(0., T_off, x, x0, para, para2);
    t_min = t_min_fun(x0[0], x0[1], para, para2);
    res_2 = v_base-v_min(0., t_min, x0[0], x0[1], para, para2);

    while ((res_1-res_2)/res_1>=0.01)
    {
        ht+=1;
        sol_1(0., T_on, x0, x, para, para2);
        sol_0(0., T_off, x, x0, para, para2);
        t_min = t_min_fun(x0[0], x0[1], para, para2);
        res_1 = res_2;
        res_2 = v_base-v_min(0., t_min, x0[0], x0[1], para, para2);
    }

    free(x0);
    free(x);

    return ht;
}

double compute_rt(double T_on, double T_off, int m, double *p0, double *para, double *para2)
{
    double t_min, res_1, res_2, v_base;
    int ht;

    double *x0 = (double*)malloc(m*sizeof(double));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }    
    double *x = (double*)malloc(m*sizeof(double));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    //
    // habituation process
    //
    for (int i = 0; i < m; i++)
    {
        x0[i] = p0[i];
    }
    ht = compute_ht(T_on, T_off, m, x0, para, para2);
    for (int i = 0; i < ht; i++)
    {
        sol_1(0., T_on, x0, x, para, para2);
        sol_0(0., T_off, x, x0, para, para2);
    }
    
    //
    // recovery time 
    // 
    double rt=0;
    double dt = 0.001;
    v_base = para[5]*(para2[1]-para2[0]);
    t_min = t_min_fun(para2[0], para2[1], para, para2);
    res_1 = v_base-v_min(0., t_min, para2[0], para2[1], para, para2);
    t_min = t_min_fun(x0[0], x0[1], para, para2);
    res_2 = v_base-v_min(0., t_min, x0[0], x0[1], para, para2);
    while ((res_1-res_2)/res_1>=0.05)
    {
        rt += dt;
        sol_0(0., rt, x0, x, para, para2);
        res_2 = v_base-v_min(0., t_min, x[0], x[1], para, para2);
    }

    free(x0);
    free(x);

    return rt;
}

int compute_potentiation(double T_on, double T_off, double *para, double *para2)
{
    double t_min, res_1, res_2, v_base, rt;
    int ht;
    int m=2;

    double *x0 = (double*)malloc(m*sizeof(double));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }
    double *x = (double*)malloc(m*sizeof(double));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    //
    // initial value
    //
    x0[0] = para2[0];
    x0[1] = para2[1];
    ht = compute_ht(T_on, T_off, m, x0, para, para2);
    rt = compute_rt(T_on, T_off, m, x0, para, para2); 
    for (int i = 0; i < ht; i++)
    {
        sol_1(0., T_on, x0, x, para, para2);
        sol_0(0., T_off, x, x0, para, para2);
    }
    sol_0(0., rt/2., x0, x, para, para2);

    //
    // second habituation process
    //
    ht = compute_ht(T_on, T_off, m, x, para, para2);
    
    free(x0);
    free(x);

    return ht;
}