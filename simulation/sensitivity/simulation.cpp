#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>

using namespace std;

# define gr (double)1.0 
# define eps (double)0.05

void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2);
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2);
void compute_crossing_time_Q(int m, int i, double *para, double *para2);
void compute_crossing_time_C1(int m, int i, double *para, double *para2);
void compute_crossing_time_C2(int m, int i, double *para, double *para2);
void plot_crossing_time_Q();
void plot_crossing_time_C1();
void plot_crossing_time_C2();
int main(int argc, const char* argv[])
{
    int n_sen = 3;
    // gf: para[0], gfs: para[1], ade: para[2], ades: para[3], Dif: para[4], Kl: para[5], rv: para[6]

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

    double para_sen[7];
    for (int i = 0; i < 7; i++)
    {
        para_sen[i] = para[i];
    }
    double noise[3] ={1.+eps, 1., 1.-eps};
    int m = 3; // Model dimension
    
    //
    //  Output the model parameters
    //
    string param_sim = "result/param/param_sim.dat";
    FILE *fparam_sim;
    fparam_sim = fopen(param_sim.c_str(), "w");
    if(fparam_sim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fprintf(fparam_sim, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n", para[0], gr, para[1], para[2], para[3], para[4], para[5], para[6]);
    fprintf(fparam_sim, "#gf\n#gr\n#gf_s\n#a_de\n#a_de_s\n#Dif\n#gam\n#rv\n");
    fclose(fparam_sim);

    string param = "result/param/param_plot.dat";
    FILE *fparam;
    fparam = fopen(param.c_str(), "w");
    if(fparam == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fprintf(fparam, "%lf\n%lf\n%lf\n%lf\n", para2[2], para2[0], para2[3], para2[1]);
    fprintf(fparam, "#Pf_s\n#Pf_u\n#Pr_s\n#Pr_u\n");
    fclose(fparam);

    string param2 = "result/param/param_plot2.dat";
    FILE *fparam2;
    fparam2 = fopen(param2.c_str(), "w");
    if(fparam2 == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fprintf(fparam2, "%lf\t%lf\n", para2[0], para2[1]);
    fprintf(fparam2, "%lf\t%lf\n", para2[2], para2[3]);
    fprintf(fparam2, "#Pf_u\t#Pr_u\n#Pf_s\t#Pr_s\n");
    fclose(fparam2);

    //
    //  Writing solutions
    //
    printf("Computing time to cross the bridge!! ");
    for (int i = 0; i < 3; i++)
    {
        para_sen[n_sen] = noise[i]*para[n_sen];
        para2[0]= ((para_sen[2]+para_sen[6]*para_sen[4])*para_sen[0]+para_sen[4]*gr)/((para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])*para_sen[2]);
        para2[1]= (para_sen[6]*para_sen[4]*para_sen[0] + (para_sen[2]+para_sen[4])*gr)/((para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])*para_sen[2]);
        para2[2]= ((para_sen[2]+para_sen[6]*para_sen[4])*para_sen[1]+para_sen[4]*gr)/((para_sen[2]+para_sen[6]*para_sen[4])*para_sen[3]+para_sen[2]*para_sen[4]);
        para2[3]= (para_sen[6]*para_sen[4]*para_sen[1] + (para_sen[3]+para_sen[4])*gr)/((para_sen[2]+para_sen[6]*para_sen[4])*para_sen[3]+para_sen[2]*para_sen[4]);
        para2[4]= -0.5*(para_sen[3]+para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])-0.5*sqrt((para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])*(para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])+4.*para_sen[6]*para_sen[4]*para_sen[4]);
        para2[5]= -0.5*(para_sen[3]+para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])+0.5*sqrt((para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])*(para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])+4.*para_sen[6]*para_sen[4]*para_sen[4]);
        para2[6]= -0.5*(para_sen[3]+para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])-0.5*sqrt((para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])*(para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])+4.*para_sen[6]*para_sen[4]*para_sen[4])+para_sen[2]+para_sen[6]*para_sen[4];
        para2[7]= para_sen[6]*para_sen[4];
        para2[8]= para_sen[4];
        para2[9]= -0.5*(para_sen[3]+para_sen[2]+para_sen[4]+para_sen[6]*para_sen[4])+0.5*sqrt((para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])*(para_sen[3]-para_sen[2]+para_sen[4]-para_sen[6]*para_sen[4])+4.*para_sen[6]*para_sen[4]*para_sen[4])+para_sen[3]+para_sen[4];
        compute_crossing_time_Q(m, i, para_sen, para2);
        compute_crossing_time_C1(m, i, para_sen, para2);
        compute_crossing_time_C2(m, i, para_sen, para2);
    }
    printf("-> Completed the computing!!\n");

    printf("Plotting result !! ");
    plot_crossing_time_Q();
    plot_crossing_time_C1();
    plot_crossing_time_C2();
    printf("-> Completed the plotting!!\n");

    return 0;
}
// analytical solution for h(t)=0
void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2)
{  
    x[0] = para2[0]+(x0[0]-x0[1]-para2[0]+para2[1])*exp(-(para[2]+para[4]+para[6]*para[4])*(t-t0))/(1.+para[6]) + (para[6]*x0[0]+x0[1]-para[6]*para2[0]-para2[1])*exp(-para[2]*(t-t0))/(1.+para[6]);
    x[1] = para2[1]-para[6]*(x0[0]-x0[1]-para2[0]+para2[1])*exp(-(para[2]+para[4]+para[6]*para[4])*(t-t0))/(1.+para[6]) + (para[6]*x0[0]+x0[1]-para[6]*para2[0]-para2[1])*exp(-para[2]*(t-t0))/(1.+para[6]);
    x[2] = x0[2] + para[5]*(para2[1]-para2[0])*(t-t0) -para[5]*(x0[0]-x0[1]-para2[0]+para2[1])*(1.-exp(-(para[2]+para[4]+para[6]*para[4])*(t-t0)))/(para[2]+para[4]+para[6]*para[4]);
}

// analytical solution for h(t)=1
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2)
{
    x[0] = para2[2] + para2[6]*(para2[9]*(x0[0]-para2[2])-para2[8]*(x0[1]-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[8]*(-para2[7]*(x0[0]-para2[2])+para2[6]*(x0[1]-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);
    x[1] = para2[3] + para2[7]*(para2[9]*(x0[0]-para2[2])-para2[8]*(x0[1]-para2[3]))*exp(para2[4]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]) + para2[9]*(-para2[7]*(x0[0]-para2[2])+para2[6]*(x0[1]-para2[3]))*exp(para2[5]*(t-t0))/(para2[6]*para2[9]-para2[8]*para2[7]);
    x[2]  = x0[2] + para[5]*(para2[3]-para2[2])*(t-t0)
        + para[5]*(para2[9]*(x0[0]-para2[2])-para2[8]*(x0[1]-para2[3]))*(para2[7]-para2[6])*(exp(para2[4]*(t-t0))-1.)/(para2[4]*(para2[6]*para2[9]-para2[8]*para2[7]))
        + para[5]*(-para2[7]*(x0[0]-para2[2])+para2[6]*(x0[1]-para2[3]))*(para2[9]-para2[8])*(exp(para2[5]*(t-t0))-1.)/(para2[5]*(para2[6]*para2[9]-para2[8]*para2[7]));
}

void compute_crossing_time_Q(int m, int i, double *para, double *para2)
{
    double T_LEARN, T_REST, t, t0; // Maximal time to compute
    double h; // time step
    int save_step, flag, day; // num. of saving step

    //
    //      Maximal time and num. of grids
    //
    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &T_LEARN);
    fscanf(fdim, "%lf", &T_REST);
    fscanf(fdim, "%lf", &h);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);
    
    double Len; // bridge length
    int TERM0, TERM1, TERM2, TERM3; // Experimental term
    string set = "setting_simu/setting_Q.dat";
    FILE *fset;
    fset = fopen(set.c_str(), "r");
    if(fset == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fset, "%lf", &Len);
    fscanf(fset, "%d", &TERM0);
    fscanf(fset, "%d", &TERM1);
    fscanf(fset, "%d", &TERM2);
    fscanf(fset, "%d", &TERM3);
    fclose(fset);

    double *x0 = (double*)malloc(m*sizeof(double *));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }
    double *x = (double*)malloc(m*sizeof(double *));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    //  iteration step
    int IT_MAX_L = T_LEARN/h;
    int IT_MAX_R = T_REST/h;
    // correcting iteration step
    if (T_LEARN >= (IT_MAX_L+1)*h)
    {
        IT_MAX_L +=1;
    }
    if (T_REST >= (IT_MAX_R+1)*h)
    {
        IT_MAX_R +=1;
    }
    
    FILE *fD;
    string cross = "result/simulation_Q/cross_time_" + std::to_string(i)+".dat";
    fD = fopen(cross.c_str(), "w");

    // substitution of initial values
    t = 0;
    day = 0;
    x[0] = para2[0];
    x[1] = para2[1];
    x[2] = 0.;

    // Experimental Phase 1 (b-Q)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
                
            }
                        
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }
    }
    
    // Experimental Phase 2 (b-N)
    for (int i = 0; i < TERM2; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
            }      
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }

    }
    
    // Experimental Phase 3 (b-Q)
    for (int i = 0; i < TERM3; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
                
            }
                        
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }

    }
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void compute_crossing_time_C1(int m, int i, double *para, double *para2)
{
    double T_LEARN, T_REST, t, t0; // Maximal time to compute
    double h; // time step
    int save_step, flag, day; // num. of saving step

    //
    //      Maximal time and num. of grids
    //
    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &T_LEARN);
    fscanf(fdim, "%lf", &T_REST);
    fscanf(fdim, "%lf", &h);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);
    
    double Len; // bridge length
    int TERM0, TERM1, TERM2, TERM3, TERM4; // Experimental term
    string set = "setting_simu/setting_C1.dat";
    FILE *fset;
    fset = fopen(set.c_str(), "r");
    if(fset == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fset, "%lf", &Len);
    fscanf(fset, "%d", &TERM0);
    fscanf(fset, "%d", &TERM1);
    fscanf(fset, "%d", &TERM2);
    fscanf(fset, "%d", &TERM3);
    fscanf(fset, "%d", &TERM4);
    fclose(fset);

    double *x0 = (double*)malloc(m*sizeof(double *));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }
    double *x = (double*)malloc(m*sizeof(double *));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    //  iteration step
    int IT_MAX_L = T_LEARN/h;
    int IT_MAX_R = T_REST/h;
    // correcting iteration step
    if (T_LEARN >= (IT_MAX_L+1)*h)
    {
        IT_MAX_L +=1;
    }
    if (T_REST >= (IT_MAX_R+1)*h)
    {
        IT_MAX_R +=1;
    }
    
    FILE *fD;
    string cross = "result/simulation_C1/cross_time_" + std::to_string(i)+".dat";
    fD = fopen(cross.c_str(), "w");

    // substitution of initial values
    t = 0.;
    day = 0;
    x[0] = para2[0];
    x[1] = para2[1];
    x[2] = 0.;

    // Experimental Phase 1 (b-N)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
                
            }
                        
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }
    }
    
    // Experimental Phase 2 (b-Q)
    for (int i = 0; i < TERM2; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
            }      
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }

    }
    
    // Experimental Phase 3 (b-N)
    for (int i = 0; i < TERM3; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
                
            }
                        
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }

    }

    // Experimental Phase 4 (b-Q)
    for (int i = 0; i < TERM4; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
            }      
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }

    }
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void compute_crossing_time_C2(int m, int i, double *para, double *para2)
{
    double T_LEARN, T_REST, t, t0; // Maximal time to compute
    double h; // time step
    int save_step, flag, day; // num. of saving step

    //
    //      Maximal time and num. of grids
    //
    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &T_LEARN);
    fscanf(fdim, "%lf", &T_REST);
    fscanf(fdim, "%lf", &h);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);
    
    double Len; // bridge length
    int TERM0, TERM1; // Experimental term
    string set = "setting_simu/setting_C2.dat";
    FILE *fset;
    fset = fopen(set.c_str(), "r");
    if(fset == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fset, "%lf", &Len);
    fscanf(fset, "%d", &TERM0);
    fscanf(fset, "%d", &TERM1);
    fclose(fset);

    double *x0 = (double*)malloc(m*sizeof(double *));
    if(x0 == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }
    double *x = (double*)malloc(m*sizeof(double *));
    if(x == NULL)
    {
        fprintf(stderr, " Failed allocating memory. \n");
        exit(1);
    }

    //  iteration step
    int IT_MAX_L = T_LEARN/h;
    int IT_MAX_R = T_REST/h;
    // correcting iteration step
    if (T_LEARN >= (IT_MAX_L+1)*h)
    {
        IT_MAX_L +=1;
    }
    if (T_REST >= (IT_MAX_R+1)*h)
    {
        IT_MAX_R +=1;
    }
    
    FILE *fD;
    string cross = "result/simulation_C2/cross_time_" + std::to_string(i)+".dat";
    fD = fopen(cross.c_str(), "w");

    // substitution of initial values
    t = 0.;
    day = 0;
    x[0] = para2[0];
    x[1] = para2[1];
    x[2] = 0.;

    // Experimental Phase 1 (b-N)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    flag =1;
                }
                
            }
                        
        }
    
        // compute resting phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = x[2];
        for(int it = 1; it <= IT_MAX_R; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            t = t0 +(double)it*h;
            
            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0+T_LEARN, x[2]);
                    flag =1;
                }   
            }             
        }
    }
    
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void plot_crossing_time_Q()
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";
    double lw_sol=5.0, lw_sol_emph = 7, lw_vel=5.0;
    int Day=9;
    
    string cross_normal = "result/simulation_Q/cross_time_1.dat";
    string cross_peps = "result/simulation_Q/cross_time_0.dat";
    string cross_meps = "result/simulation_Q/cross_time_2.dat";

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    
    fprintf(gpid, "reset\n");
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_Q_sen.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", 0.5, 13.);
    fprintf(gpid, "set xtics %d, %d, %d nomirror scale 1.5\n",0, 1, Day);
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set tics font '%s''%s' out\n", font_tics,",25");
    fprintf(gpid, "set object rect from 0.5,0.5 to 6.5,13 fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n");
    fprintf(gpid, "set object rect from 6.5,0.5 to 8.5,13 fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n");
    fprintf(gpid, "set object rect from 8.5,0.5 to 9.5,13 fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n");

    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2 w l lw 4  lc rgb 'black', ",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 2 lc rgb 'black',",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4 dt (6,3) lc rgb 'red', ",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 8 lc rgb 'red',",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4 dt (2,2) lc rgb 'blue', ",cross_meps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 10 lc rgb 'blue'\n",cross_meps.c_str());

    fflush(gpid);
    pclose(gpid);
}

void plot_crossing_time_C1()
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";
    double lw_sol=5.0, lw_sol_emph = 7, lw_vel=5.0;
    int Day=9;
    
    string cross_normal = "result/simulation_C1/cross_time_1.dat";
    string cross_peps = "result/simulation_C1/cross_time_0.dat";
    string cross_meps = "result/simulation_C1/cross_time_2.dat";

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    
    fprintf(gpid, "reset\n");
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_C1_sen.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", 0.5, 13.);
    fprintf(gpid, "set xtics %d, %d, %d nomirror scale 1.5\n",0, 1, Day);
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set tics font '%s''%s' out\n", font_tics,",25");
    fprintf(gpid, "set object rect from 0.5,0.5 to 5.5,13 fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n");
    fprintf(gpid, "set object rect from 5.5,0.5 to 6.5,13 fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n");
    fprintf(gpid, "set object rect from 6.5,0.5 to 8.5,13 fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n");
    fprintf(gpid, "set object rect from 8.5,0.5 to 9.5,13 fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n");

    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2 w l lw 4  lc rgb 'black', ",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 2 lc rgb 'black',",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4 dt (6,3) lc rgb 'red', ",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 8 lc rgb 'red',",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4 dt (2,2) lc rgb 'blue', ",cross_meps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3 pt 10 lc rgb 'blue'\n",cross_meps.c_str());

    fflush(gpid);
    pclose(gpid);
}

void plot_crossing_time_C2()
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";
    int Day=9;
    
    string cross_normal = "result/simulation_C2/cross_time_1.dat";
    string cross_peps = "result/simulation_C2/cross_time_0.dat";
    string cross_meps = "result/simulation_C2/cross_time_2.dat";

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    
    fprintf(gpid, "reset\n");
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_C2_sen.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", 0.5, 13.);
    fprintf(gpid, "set xtics %d, %d, %d nomirror scale 1.5\n",0, 1, Day);
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set tics font '%s''%s' out\n", font_tics,",25");
    fprintf(gpid, "set object rect from 0.5,0.5 to 9.5,13 fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n");

    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2 w l lw 4.  lc rgb 'black', ",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3. pt 2 lc rgb 'black',",cross_normal.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4. dt (6,3) lc rgb 'red', ",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3. pt 8 lc rgb 'red',",cross_peps.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 4. dt (2,2) lc rgb 'blue', ",cross_meps.c_str());
    fprintf(gpid, "'%s' using 1:2 w p ps 3. lw 3. pt 10 lc rgb 'blue'\n",cross_meps.c_str());

    fflush(gpid);
    pclose(gpid);
}
