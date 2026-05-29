#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>

using namespace std;

# define gr (double)1.0 

void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2);
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2);
void compute_crossing_time_Q(int m, double *para, double *para2);
void compute_crossing_time_C1(int m, double *para, double *para2);
void compute_crossing_time_C2(int m, double *para, double *para2);
void plot_crossing_time_Q(double y_min, double y_max);
void plot_crossing_time_C1(double y_min, double y_max);
void plot_crossing_time_C2(double y_min, double y_max);
void plot_time_course_Q(int day, double p_min, double p_max, double l_min, double l_max, double v_min, double v_max, double m_min, double m_max,double *para, double *para2);

int main(int argc, const char* argv[])
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
    compute_crossing_time_Q(m, para, para2);
    compute_crossing_time_C1(m, para, para2);
    compute_crossing_time_C2(m, para, para2);
    printf("-> Completed the computing!!\n");

    double y_min, y_max, p_min, p_max, l_min, l_max, v_min, v_max, m_min, m_max;
    int day_phase;
    printf("Plotting result !! ");
    y_min = 0.9;
    y_max = 4.3;
    day_phase = 8;
    p_min = para2[0]-0.2; //para2[0]=Pf_u
    p_max = para2[3];      //para2[3]=Pr_s
    l_min = -0.01;
    l_max = 0.28;
    v_min = 0.0;
    v_max = 0.012; 
    m_min = 81.5;
    m_max = 85.; 
    plot_crossing_time_Q(y_min, y_max);
    plot_crossing_time_C1(y_min, y_max);
    plot_crossing_time_C2(y_min, y_max);
    plot_time_course_Q(day_phase, p_min, p_max, l_min, l_max, v_min, v_max, m_min, m_max, para, para2);
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

void compute_crossing_time_Q(int m, double *para, double *para2)
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
    
    FILE *fw, *fD, *fi;
    string data = "result/simulation_Q/data.dat";
    string cross = "result/simulation_Q/cross_time.dat";
    string ini_con = "result/simulation_Q/ini_con.dat";
    fD = fopen(cross.c_str(), "w");
    fclose(fD);
    fi = fopen(ini_con.c_str(), "w");
    fclose(fi);

    // substitution of initial values
    t = -(T_LEARN+T_REST)*TERM0;
    day = -TERM0;
    x0[0] = para2[0];
    x0[1] = para2[1];
    x0[2] = 0.;

    // Output of the initial state
    fw = fopen(data.c_str(), "w");
    fprintf(fw, "%lf\t", t);
    for(int i = 0; i < m; i++)
    {
        fprintf(fw, "%lf\t", x[i]);              
    }
    fprintf(fw, "\n");

    // Dummy Phase (b-N)
    for (int i = 0; i < TERM0; i++)
    {
        day +=1;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        for(int it = 1; it <= IT_MAX_L+IT_MAX_R; it++)
        {   
            sol_0(t0, t, x0, x, para, para2);

            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                fprintf(fw, "%lf\t", x[0]);
                fprintf(fw, "%lf\t", x[1]);
                fprintf(fw, "%lf\t", 0.);
                fprintf(fw, "\n");
            }           
        }

    }

    // Experimental Phase 1 (b-Q)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
                     
            }            
        }
    }
    
    // Experimental Phase 2 (b-N)
    for (int i = 0; i < TERM2; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

        }

    }
    
    // Experimental Phase 3 (b-Q)
    for (int i = 0; i < TERM3; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

        }

    }
    
    fprintf(fw, "# time\tpf\tpr\tl\n");
    fclose(fw);
    fD = fopen(cross.c_str(), "a");
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void compute_crossing_time_C1(int m, double *para, double *para2)
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
    
    FILE *fw, *fD, *fi;
    string data = "result/simulation_C1/data.dat";
    string cross = "result/simulation_C1/cross_time.dat";
    string ini_con = "result/simulation_C1/ini_con.dat";
    fD = fopen(cross.c_str(), "w");
    fclose(fD);
    fi = fopen(ini_con.c_str(), "w");
    fclose(fi);

    // substitution of initial values
    t = -(T_LEARN+T_REST)*TERM0;
    day = -TERM0;

    // Output of the initial state
    fw = fopen(data.c_str(), "w");
    fprintf(fw, "%lf\t", t);
    for(int i = 0; i < m; i++)
    {
        fprintf(fw, "%lf\t", x[i]);              
    }
    fprintf(fw, "\n");

    // Dummy Phase (b-N)
    for (int i = 0; i < TERM0; i++)
    {
        day +=1;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = para2[0];
        x0[1] = para2[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L+IT_MAX_R; it++)
        {   
            sol_0(t0, t, x0, x, para, para2);

            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }           
        }

    }

    // Experimental Phase 1 (b-N)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
                     
            }            
        }
    }
    
    // Experimental Phase 2 (b-Q)
    for (int i = 0; i < TERM2; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

        }

    }
    
    // Experimental Phase 3 (b-N)
    for (int i = 0; i < TERM3; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

        }

    }
    
     // Experimental Phase 4 (b-Q)
    for (int i = 0; i < TERM4; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_1(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

        }

    }
    
    fprintf(fw, "# time\tpf\tpr\tl\n");
    fclose(fw);
    fD = fopen(cross.c_str(), "a");
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void compute_crossing_time_C2(int m, double *para, double *para2)
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
    
    FILE *fw, *fD, *fi;
    string data = "result/simulation_C2/data.dat";
    string cross = "result/simulation_C2/cross_time.dat";
    string ini_con = "result/simulation_C2/ini_con.dat";
    fD = fopen(cross.c_str(), "w");
    fclose(fD);
    fi = fopen(ini_con.c_str(), "w");
    fclose(fi);

    // substitution of initial values
    t = -(T_LEARN+T_REST)*TERM0;
    day = -TERM0;

    // Output of the initial state
    fw = fopen(data.c_str(), "w");
    fprintf(fw, "%lf\t", t);
    for(int i = 0; i < m; i++)
    {
        fprintf(fw, "%lf\t", x[i]);              
    }
    fprintf(fw, "\n");

    // Dummy Phase (b-N)
    for (int i = 0; i < TERM0; i++)
    {
        day +=1;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = para2[0];
        x0[1] = para2[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L+IT_MAX_R; it++)
        {   
            sol_0(t0, t, x0, x, para, para2);

            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }           
        }

    }

    // Experimental Phase 1 (b-N)
    for (int i = 0; i < TERM1; i++)
    {
        day +=1;
        flag =0;

        // Output of the initial state
        fclose(fw);
        fi = fopen(ini_con.c_str(), "a");
        fprintf(fi, "%d\t%lf\t%lf\t", day, x0[0], x0[1]);
        fprintf(fi, "\n");
        fclose(fi);
        fw = fopen(data.c_str(), "a");

        // compute learning phase
        t0 = t;
        x0[0] = x[0];
        x0[1] = x[1];
        x0[2] = 0.;
        for(int it = 1; it <= IT_MAX_L; it++)
        {
            sol_0(t0, t, x0, x, para, para2);
            
            //  Update
            t = t0 +(double)it*h;
            
            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
            }

            if (flag == 0)
            {
                if (x[2] >= Len)
                {
                    fclose(fw);
                    fD = fopen(cross.c_str(), "a");
                    fprintf(fD, "%d\t%lf\t%lf\n", day, t-t0, x[2]);
                    fclose(fD);
                    fw = fopen(data.c_str(), "a");
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
            
            //  Update
            t = t0 +(double)it*h;

            if(it % save_step == 0) // Output the solution at each save_step
            {
                fprintf(fw, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw, "%lf\t", x[i]);              
                }
                fprintf(fw, "\n");
                     
            }            
        }
    }
    
    fprintf(fw, "# time\tpf\tpr\tl\n");
    fclose(fw);
    fD = fopen(cross.c_str(), "a");
    fprintf(fD, "# day\ttime\tl\n");
    fclose(fD);

    free(x);     
    free(x0);
}

void plot_crossing_time_Q(double y_min, double y_max)
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";

    string cross_mean_CI_data = "data/Q_mean_CI_data.dat";
    string cross_simu = "result/simulation_Q/cross_time.dat";

    int Day = 9;

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // Open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_Q.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    //fprintf(gpid, "set xlabel '%s'\n", "day");
    //fprintf(gpid, "set ylabel '%s'\n", "crossing time [h]");
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set xtics %d, %d, %d\n",0, 1, Day);
    fprintf(gpid, "set tics font '%s''%s'\n", font_tics,",25");
    
    fprintf(gpid, "set object rect from 0.5,%lf to 6.5,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", y_min, y_max); 
    fprintf(gpid, "set object rect from 6.5,%lf to 8.5,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", y_min, y_max);
    fprintf(gpid, "set object rect from 8.5,%lf to 9.5,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", y_min, y_max);
    
    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2:3:4 with yerrorbars pt 0 lw 3. lc rgb 'black',",cross_mean_CI_data.c_str()); // plot error bar
    fprintf(gpid, "'%s' using 1:2 with points ps 2. lw 2.5 pt 6 lc rgb 'black', ",cross_mean_CI_data.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 2.5 lc rgb '#FF4040', ",cross_simu.c_str()); // plot line
    fprintf(gpid, "'%s' using 1:2 w p ps 3.5 lw 5. pt 2 lc rgb 'red'\n",cross_simu.c_str()); // plot point
    fflush(gpid);
    pclose(gpid);
}

void plot_crossing_time_C1(double y_min, double y_max)
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";

    string cross_mean_CI_data = "data/C1_mean_CI_data.dat";
    string cross_simu = "result/simulation_C1/cross_time.dat";

    int Day = 9;

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // Open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_C1.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    //fprintf(gpid, "set xlabel '%s'\n", "day");
    //fprintf(gpid, "set ylabel '%s'\n", "crossing time [h]");
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set xtics %d, %d, %d\n",0, 1, Day);
    fprintf(gpid, "set tics font '%s''%s'\n", font_tics,",25");
    
    fprintf(gpid, "set object rect from 0.5,%lf to 5.5,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", y_min, y_max); 
    fprintf(gpid, "set object rect from 5.5,%lf to 6.5,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", y_min, y_max);
    fprintf(gpid, "set object rect from 6.5,%lf to 8.5,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", y_min, y_max);
    fprintf(gpid, "set object rect from 8.5,%lf to 9.5,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", y_min, y_max);

    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2:3:4 with yerrorbars pt 0 lw 3. lc rgb 'black',",cross_mean_CI_data.c_str()); // plot error bar
    fprintf(gpid, "'%s' using 1:2 with points ps 2. lw 2.5 pt 6 lc rgb 'black', ",cross_mean_CI_data.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 2.5 lc rgb '#FF4040', ",cross_simu.c_str()); // plot line
    fprintf(gpid, "'%s' using 1:2 w p ps 3.5 lw 5. pt 2 lc rgb 'red'\n",cross_simu.c_str()); // plot point
    fflush(gpid);
    pclose(gpid);
}

void plot_crossing_time_C2(double y_min, double y_max)
{
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";

    string cross_mean_CI_data = "data/C2_mean_CI_data.dat";
    string cross_simu = "result/simulation_C2/cross_time.dat";

    int Day = 9;

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // Open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/cross_time_C2.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", 1.);
    //fprintf(gpid, "set xlabel '%s'\n", "day");
    //fprintf(gpid, "set ylabel '%s'\n", "crossing time [h]");
    fprintf(gpid, "set xrange [%f:%f]\n", 0.5, Day+0.5);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set xtics %d, %d, %d\n",0, 1, Day);
    fprintf(gpid, "set tics font '%s''%s'\n", font_tics,",25");
    
    fprintf(gpid, "set object rect from 0.5,%lf to 9.5,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", y_min, y_max); 

    fprintf(gpid,"set border front lw 2\n");
    fprintf(gpid, "plot '%s' using 1:2:3:4 with yerrorbars pt 0 lw 3. lc rgb 'black',",cross_mean_CI_data.c_str()); // plot error bar
    fprintf(gpid, "'%s' using 1:2 with points ps 2. lw 2.5 pt 6 lc rgb 'black', ",cross_mean_CI_data.c_str());
    fprintf(gpid, "'%s' using 1:2 w l lw 2.5 lc rgb '#FF4040', ",cross_simu.c_str()); // plot line
    fprintf(gpid, "'%s' using 1:2 w p ps 3.5 lw 5. pt 2 lc rgb 'red'\n",cross_simu.c_str()); // plot point
    fflush(gpid);
    pclose(gpid);
}

void plot_time_course_Q(int day, double p_min, double p_max, double l_min, double l_max, double v_min, double v_max, double m_min, double m_max, double *para, double *para2)
{
    double lw_sol=4.0, lw_sol_emph = 6.0, lw_vel=5.0, lw_acc=4.0, lw_time=6.0;
    int N_acc = 4;
    char font_title[] = "Helvetica,20", font_label[] = "Helvetica,20", font_tics[] = "Helvetica";

    double T_LEARN, T_REST, IT_STEP;
    int save_step; 

    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &T_LEARN);
    fscanf(fdim, "%lf", &T_REST);
    fscanf(fdim, "%lf", &IT_STEP);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);

    double Len; // bridge length
    int TERM0, TERM1, TERM2, TERM3, Day; // Experimental term
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
    Day = TERM1+TERM2+TERM3;

    double Pf_u=para2[0], Pr_u=para2[1], Pf_s=para2[2], Pr_s=para2[3];
    double x_min = Pf_u-0.1, x_max = Pf_s+0.1;
    double y_min = Pr_s-(Pf_s-Pf_u)+0.1, y_max = Pr_s+0.3;
    int step1 = (T_LEARN+T_REST)/(IT_STEP*save_step);

    string data = "result/simulation_Q/data.dat";
    string equi = "result/param/param_plot2.dat";

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // Open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }

    //
    // plot phase diagram 
    //
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/pf_pr_Q.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", (y_max-y_min)/(x_max-x_min));
    fprintf(gpid, "set xrange [%f:%f]\n", x_min, x_max);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set xtics nomirror scale 1.5\n");
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set tics font '%s'\n", font_tics);
    fprintf(gpid,"set border back lw 2\n");
    
    // plot iso-velocity lines
    fprintf(gpid, "plot [%f:%f] x w l lw %f lc rgb 'red' dt (16,8),",x_min, x_max, lw_vel);
    fprintf(gpid, "[%f:%f] x+%f w l lw %f lc rgb 'red' dt (16,8),",x_min, x_max, Pr_u-Pf_u,lw_vel);  

    // plot iso-acceleration lines
    double k = (para[3]+para[4]+para[6]*para[2])/(para[2]+para[4]+para[6]*para[2]);
    for (int i = 0; i < N_acc; i++)
    {
        fprintf(gpid, "[%f:%f] %lf*(x-%lf)+%lf w l lw %f lc rgb '%s' dt (4,4) notitle",x_min, x_max, k, Pf_u+i*0.4*(Pf_s-Pf_u)/(N_acc-1), Pr_u, lw_acc, "blue");
        if (i!=N_acc-1)
        {
            fprintf(gpid,",");
        }
    }
    fprintf(gpid,",");

    // plot solution
    if (day < 7)
    {
        fprintf(gpid, "'%s' using 2:3 every ::%d::%d w l lw %f lc rgb 'dark-grey',",data.c_str(), 1*step1, (day+1)*step1, lw_sol); // plot day1-day6 (gray)
    }
    else
    {
        fprintf(gpid, "'%s' using 2:3 every ::%d::%d w l lw %f lc rgb 'dark-grey',",data.c_str(), 1*step1, 7*step1, lw_sol); // plot day1-day6 (gray)
        fprintf(gpid, "'%s' using 2:3 every ::%d::%d w l lw %f lc rgb 'black',",data.c_str(), 7*step1, (day+1)*step1, lw_sol_emph);// plot day7-day9 (black)

    }
    
    // plot equilibrium point
    fprintf(gpid, "'%s' using 1:2 every ::1::1 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str()); 
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str()); 
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 3.2 lc rgb 'white'\n",equi.c_str());

    //
    // plot time course of pf, pr
    //
    int n;
    fprintf(gpid, "reset\n"); 
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/time_p_Q.png'\n");
    fprintf(gpid, "set xrange [%f:%f]\n", -(T_LEARN+T_REST)/2., (T_LEARN+T_REST)*Day);
    fprintf(gpid, "set yrange [%f:%f]\n", p_min, p_max);
    fprintf(gpid,"set border lw 2\n");
    fprintf(gpid, "set xtics (");
    for (int i = 0; i <= Day; i++)
    {   
        fprintf(gpid, "%f", (T_LEARN+T_REST)*i);
        if (i!=Day)
        {
            fprintf(gpid, ",");
        }
    }
    fprintf(gpid, ") font '%s''%s'\n", font_tics, ",20");
    for (int i = 0; i < TERM0; i++)
    {
        n = i-TERM0;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, p_min, (T_LEARN+T_REST)*(n+1), p_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), p_min, (T_LEARN+T_REST)*(n+1), p_max, 5.);
    }
    for (int i = 0; i < TERM1; i++)
    {
        n = i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, p_min, T_LEARN+(T_LEARN+T_REST)*n, p_max); 
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, p_min, (T_LEARN+T_REST)*(n+1), p_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), p_min, (T_LEARN+T_REST)*(n+1), p_max, 5.);
 
    }
    for (int i = 0; i < TERM2; i++)
    {
        n = TERM1+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, p_min, (T_LEARN+T_REST)*(n+1), p_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), p_min, (T_LEARN+T_REST)*(n+1), p_max, 5.);
 
    }
    for (int i = 0; i < TERM3; i++)
    {
        n = TERM1+TERM2+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, p_min, T_LEARN+(T_LEARN+T_REST)*n, p_max); 
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, p_min, (T_LEARN+T_REST)*(n+1), p_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), p_min, (T_LEARN+T_REST)*(n+1), p_max, 5.);
 
    }
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "plot '%s' using 1:2 w l lw %lf lc rgb 'red' title 'p_{f}',",data.c_str(), lw_time);
    fprintf(gpid, "'%s' using 1:3 w l lw %lf lc rgb 'blue' title 'p_{b}'\n",data.c_str(), lw_time);

    //
    // plot time course of l
    //
    fprintf(gpid, "reset\n"); 
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/time_l_Q.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set xrange [%f:%f]\n", -(T_LEARN+T_REST)/2., (T_LEARN+T_REST)*Day);
    fprintf(gpid, "set yrange [%f:%f]\n", l_min, l_max);
    fprintf(gpid, "set xtics (");
    for (int i = 0; i <= Day; i++)
    {
        fprintf(gpid, "%f", (T_LEARN+T_REST)*i);
        if (i!=Day)
        {
            fprintf(gpid, ",");
        }
    }
    fprintf(gpid, ") font '%s''%s'\n", font_tics,",20");

    for (int i = 0; i < TERM0; i++)
    {
        n = i-TERM0;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, l_min, (T_LEARN+T_REST)*(n+1), l_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), l_min, (T_LEARN+T_REST)*(n+1), l_max, 5.);

    }
    for (int i = 0; i < TERM1; i++)
    {
        n = i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, l_min, T_LEARN+(T_LEARN+T_REST)*n, l_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, l_min, (T_LEARN+T_REST)*(n+1), l_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), l_min, (T_LEARN+T_REST)*(n+1), l_max, 5.);
    }
    for (int i = 0; i < TERM2; i++)
    {
        n = TERM1+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, l_min, (T_LEARN+T_REST)*(n+1), l_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), l_min, (T_LEARN+T_REST)*(n+1), l_max, 5.);
    }
    for (int i = 0; i < TERM3; i++)
    {
        n = TERM1+TERM2+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, l_min, T_LEARN+(T_LEARN+T_REST)*n, l_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, l_min, (T_LEARN+T_REST)*(n+1), l_max);
        fprintf(gpid,"set arrow from %lf,%lf to %lf,%lf front nohead dt (5,2.5) linewidth %lf\n", (T_LEARN+T_REST)*(n+1), l_min, (T_LEARN+T_REST)*(n+1), l_max, 5.);
    }
    fprintf(gpid, "set border lw 2\n");
    fprintf(gpid, "plot '%s' using 1:4 w l lw %lf lc 1\n",data.c_str(), lw_time);

    //
    // plot time course of velocity
    //
    fprintf(gpid, "reset\n"); 
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/time_velocity_Q.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set xrange [%f:%f]\n", -(T_LEARN+T_REST)/2., (T_LEARN+T_REST)*Day);
    fprintf(gpid, "set yrange [%f:%f]\n", v_min, v_max);
    fprintf(gpid, "set xtics (");
    for (int i = 0; i <= Day; i++)
    {
        fprintf(gpid, "%f", (T_LEARN+T_REST)*i);
        if (i!=Day)
        {
            fprintf(gpid, ",");
        }
    }
    fprintf(gpid, ") font '%s''%s'\n", font_tics,",20");

    for (int i = 0; i < TERM0; i++)
    {
        n = i-TERM0;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, v_min, (T_LEARN+T_REST)*(n+1), v_max);
    }
    for (int i = 0; i < TERM1; i++)
    {
        n = i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, v_min, T_LEARN+(T_LEARN+T_REST)*n, v_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, v_min, (T_LEARN+T_REST)*(n+1), v_max);
    }
    for (int i = 0; i < TERM2; i++)
    {
        n = TERM1+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, v_min, (T_LEARN+T_REST)*(n+1), v_max);
    }
    for (int i = 0; i < TERM3; i++)
    {
        n = TERM1+TERM2+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, v_min, T_LEARN+(T_LEARN+T_REST)*n, v_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, v_min, (T_LEARN+T_REST)*(n+1), v_max);
    }
    fprintf(gpid, "set border lw 2\n");
    fprintf(gpid, "plot '%s' using 1:(%lf*($3-$2)) w l lw %lf lc 1\n",data.c_str(), para[5], lw_time);

    //
    // plot time course of m=pf+pr
    //
    fprintf(gpid, "reset\n"); 
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/time_m_Q.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set xrange [%f:%f]\n", -(T_LEARN+T_REST)/2., (T_LEARN+T_REST)*Day);
    fprintf(gpid, "set yrange [%f:%f]\n", m_min, m_max);
    fprintf(gpid, "set xtics (");
    for (int i = 0; i <= Day; i++)
    {
        fprintf(gpid, "%f", (T_LEARN+T_REST)*i);
        if (i!=Day)
        {
            fprintf(gpid, ",");
        }
    }
    fprintf(gpid, ") font '%s''%s'\n", font_tics,",20");

    for (int i = 0; i < TERM0; i++)
    {
        n = i-TERM0;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, m_min, (T_LEARN+T_REST)*(n+1), m_max);
    }
    for (int i = 0; i < TERM1; i++)
    {
        n = i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, m_min, T_LEARN+(T_LEARN+T_REST)*n, m_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, m_min, (T_LEARN+T_REST)*(n+1), m_max);
    }
    for (int i = 0; i < TERM2; i++)
    {
        n = TERM1+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, m_min, (T_LEARN+T_REST)*(n+1), m_max);
    }
    for (int i = 0; i < TERM3; i++)
    {
        n = TERM1+TERM2+i;
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'black' fillstyle solid 0.2 noborder behind\n", (T_LEARN+T_REST)*n, m_min, T_LEARN+(T_LEARN+T_REST)*n, m_max);
        fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", T_LEARN+(T_LEARN+T_REST)*n, m_min, (T_LEARN+T_REST)*(n+1), m_max);
    }
    fprintf(gpid, "set border lw 2\n");
    fprintf(gpid, "plot '%s' using 1:($2+$3) w l lw %lf lc 1\n",data.c_str(), lw_time);

    fflush(gpid);
    pclose(gpid);
}