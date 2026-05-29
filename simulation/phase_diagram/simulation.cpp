#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h>

using namespace std;

# define gr (double)1.0 

void sol_0(double t0, double t, double *x0, double *x, double *para, double *para2);
void sol_1(double t0, double t, double *x0, double *x, double *para, double *para2);
void compute_0(int m, double x_min, double x_max, double y_min, double y_max, double *para, double *para2);
void compute_1(int m, double x_min, double x_max, double y_min, double y_max, double *para, double *para2);
void plot(double x_min, double x_max, double y_min, double y_max, double *para, double *para2);


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

    int m = 2; // dimension

    //
    //  Output the model parameters
    //
    string param2 = "result/param/param_plot.dat";
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
    double Pf_u=para2[0], Pr_u=para2[1], Pf_s=para2[2], Pr_s=para2[3];
    double x_min = Pf_u-0.1, x_max = Pf_s+0.1;
    double y_min = Pr_s-(Pf_s-Pf_u)+0.1, y_max = Pr_s+0.3;

    compute_0(m, x_min, x_max, y_min, y_max, para, para2);
    compute_1(m, x_min, x_max, y_min, y_max, para, para2);
    plot(x_min, x_max, y_min, y_max, para, para2);
    
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

void compute_0(int m, double x_min, double x_max, double y_min, double y_max, double *para, double *para2)
{
    double tf, IT_STEP;
    int save_step;
    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &tf);
    fscanf(fdim, "%lf", &IT_STEP);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);

    double x_h_ini, y_h_ini;
    int x_d_ini, y_d_ini;
    string ini = "setting_simu/initial.dat";
    FILE *fini;
    fini = fopen(ini.c_str(), "r");
    if(fini == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fini, "%d", &x_d_ini);
    fscanf(fini, "%d", &y_d_ini);
    fclose(fini);

    //  iteration step
    int IT_MAX = tf/IT_STEP;
    x_h_ini = (x_max-x_min)/((double)x_d_ini);
    y_h_ini = (y_max-y_min)/((double)y_d_ini);
    // correcting iteration step
    if (tf >= (IT_MAX+1)*IT_STEP)
    {
        IT_MAX +=1;
    }

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
    
    int data_num = 0;
    double t;

    FILE *fw0;
    string data0 = "result/data/data0.dat";
    fw0 = fopen(data0.c_str(), "w");
    
    // pf_0 = some values, pr_0 = y_min (fixed)
    for (int i = 0; i <= x_d_ini; i++)
    {
        t = 0.;

        // substitution of initailal values
        x0[0] = x_min+(double)(i*x_h_ini);
        x0[1] = y_min;

        // output of initailal values
        fprintf(fw0, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw0, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw0, "\t");                
        }
        fprintf(fw0, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_0(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;

            if(it % save_step == 0)
            {
                fprintf(fw0, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw0, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw0, "\t");                
                }
                fprintf(fw0, "\n");        
            }        
        }
        fprintf(fw0, "\n");

        data_num +=1;
    }
    
    // pf_0 = some values, pr_0 = y_max (fixed)
    for (int i = 0; i <= x_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_min+(double)(i*x_h_ini);
        x0[1] = y_max;

        // output of initailal values
        fprintf(fw0, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw0, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw0, "\t");                
        }
        fprintf(fw0, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_0(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw0, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw0, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw0, "\t");                
                }
                fprintf(fw0, "\n");        
            }        
        }
        fprintf(fw0, "\n");

        data_num +=1;
    }
    
    // pf_0 = x_min (fixed), pr_0 = some values
    for (int i = 1; i < y_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_min;
        x0[1] = y_min+(double)(i*y_h_ini);

        // output of initailal values
        fprintf(fw0, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw0, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw0, "\t");                
        }
        fprintf(fw0, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_0(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw0, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw0, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw0, "\t");                
                }
                fprintf(fw0, "\n");        
            }        
        }
        fprintf(fw0, "\n");

        data_num +=1;
    }
    
    // pf_0 = x_max (fixed), pr_0 = some values
    for (int i = 1; i < y_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_max;
        x0[1] = y_min+(double)(i*y_h_ini);

        // output of initailal values
        fprintf(fw0, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw0, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw0, "\t");                
        }
        fprintf(fw0, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_0(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw0, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw0, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw0, "\t");                
                }
                fprintf(fw0, "\n");        
            }        
        }
        fprintf(fw0, "\n");

        data_num +=1;
    }

    fprintf(fw0, "#time\tpf\tpr\n");
    fclose(fw0);

    string param = "result/data/param_plot0.dat";
    FILE *fparam;
    fparam = fopen(param.c_str(), "w");
    if(fparam == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    
    fprintf(fparam, "%d\n", data_num);
    fprintf(fparam, "#data_num\n");
    fclose(fparam);

    free(x0);
    free(x);
}

void compute_1(int m, double x_min, double x_max, double y_min, double y_max, double *para, double *para2)
{
    double tf, IT_STEP;
    int save_step;
    string dim = "setting_simu/dim.dat";
    FILE *fdim;
    fdim = fopen(dim.c_str(), "r");
    if(fdim == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fdim, "%lf", &tf);
    fscanf(fdim, "%lf", &IT_STEP);
    fscanf(fdim, "%d", &save_step);
    fclose(fdim);

    double x_h_ini, y_h_ini;
    int x_d_ini, y_d_ini;
    string ini = "setting_simu/initial.dat";
    FILE *fini;
    fini = fopen(ini.c_str(), "r");
    if(fini == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fini, "%d", &x_d_ini);
    fscanf(fini, "%d", &y_d_ini);
    fclose(fini);

    //  iteration step
    int IT_MAX = tf/IT_STEP;
    x_h_ini = (x_max-x_min)/((double)x_d_ini);
    y_h_ini = (y_max-y_min)/((double)y_d_ini);
    // correcting iteration step
    if (tf >= (IT_MAX+1)*IT_STEP)
    {
        IT_MAX +=1;
    }

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
    
    int data_num = 0;
    double t;

    FILE *fw1;
    string data1 = "result/data/data1.dat";
    fw1 = fopen(data1.c_str(), "w");
    
    // pf_0 = some values, pr_0 = y_min (fixed)
    for (int i = 0; i <= x_d_ini; i++)
    {
        t = 0.;

        // substitution of initailal values
        x0[0] = x_min+(double)(i*x_h_ini);
        x0[1] = y_min;

        // output of initailal values
        fprintf(fw1, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw1, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw1, "\t");                
        }
        fprintf(fw1, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_1(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;

            if(it % save_step == 0)
            {
                fprintf(fw1, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw1, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw1, "\t");                
                }
                fprintf(fw1, "\n");        
            }        
        }
        fprintf(fw1, "\n");

        data_num +=1;
    }
    
    // pf_0 = some values, pr_0 = y_max (fixed)
    for (int i = 0; i <= x_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_min+(double)(i*x_h_ini);
        x0[1] = y_max;

        // output of initailal values
        fprintf(fw1, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw1, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw1, "\t");                
        }
        fprintf(fw1, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_1(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw1, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw1, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw1, "\t");                
                }
                fprintf(fw1, "\n");        
            }        
        }
        fprintf(fw1, "\n");

        data_num +=1;
    }
    
    // pf_0 = x_min (fixed), pr_0 = some values
    for (int i = 1; i < y_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_min;
        x0[1] = y_min+(double)(i*y_h_ini);

        // output of initailal values
        fprintf(fw1, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw1, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw1, "\t");                
        }
        fprintf(fw1, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_1(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw1, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw1, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw1, "\t");                
                }
                fprintf(fw1, "\n");        
            }        
        }
        fprintf(fw1, "\n");

        data_num +=1;
    }
    
    // pf_0 = x_max (fixed), pr_0 = some values
    for (int i = 1; i < y_d_ini; i++)
    {
        t = 0.;
        // substitution of initailal values
        x0[0] = x_max;
        x0[1] = y_min+(double)(i*y_h_ini);

        // output of initailal values
        fprintf(fw1, "%lf\t", t);
        for(int i = 0; i < m; i++)
        {
            fprintf(fw1, "%lf", x0[i]);
            if(i < m-1)
                fprintf(fw1, "\t");                
        }
        fprintf(fw1, "\n");
        
        // calculation and result writing
        for(int it = 1; it <= IT_MAX; it++)
        {
            sol_1(0., t, x0, x, para, para2);
            t = (double)it*IT_STEP;
            
            if(it % save_step == 0)
            {
                fprintf(fw1, "%lf\t", t);
                for(int i = 0; i < m; i++)
                {
                    fprintf(fw1, "%lf", x[i]);
                    if(i < m-1)
                        fprintf(fw1, "\t");                
                }
                fprintf(fw1, "\n");        
            }        
        }
        fprintf(fw1, "\n");

        data_num +=1;
    }

    fprintf(fw1, "#time\tpf\tpr\n");
    fclose(fw1);

    string param = "result/data/param_plot1.dat";
    FILE *fparam;
    fparam = fopen(param.c_str(), "w");
    if(fparam == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    
    fprintf(fparam, "%d\n", data_num);
    fprintf(fparam, "#data_num\n");
    fclose(fparam);

    free(x0);
    free(x);
}

void plot(double x_min, double x_max, double y_min, double y_max, double *para, double *para2)
{
    double Pf_u=para2[0], Pr_u=para2[1], Pf_s=para2[2], Pr_s=para2[3];
    int data0_num, data1_num;

    string param0 = "result/data/param_plot0.dat";
    FILE *fparam0;
    fparam0 = fopen(param0.c_str(), "r");
    if(fparam0 == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fparam0, "%d", &data0_num);
    fclose(fparam0);

    string param1 = "result/data/param_plot1.dat";
    FILE *fparam1;
    fparam1 = fopen(param1.c_str(), "r");
    if(fparam1 == NULL) {
        printf(" Couldn't open files.\n");
        exit(1);
    }
    fscanf(fparam1, "%d", &data1_num);
    fclose(fparam1);

    string data0 = "result/data/data0.dat";
    string data1 = "result/data/data1.dat";
    string equi = "result/param/param_plot.dat";
    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // open gnuplot
    if(gpid == NULL)
    {
        printf("I can't find gnuplot.\n"); exit(1);
    }
    
    double lw_sol=4.0, lw_vel=5.0;
    char font_title[] = "Arial,20", font_label[] = "Arial,20", font_tics[] = "Arial";

    // h(t)=1
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/Phase1.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", (y_max-y_min)/(x_max-x_min));
    fprintf(gpid, "set xrange [%f:%f]\n", x_min, x_max);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set tics font '%s'\n", font_tics);
    fprintf(gpid,"set border back lw 2\n");
    fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", x_min, y_min, x_max, y_max);
    fprintf(gpid, "plot '%s' using 2:3 every :::%d::%d w l lw %f lc rgb 'black', " , data1.c_str(), 0, 0, lw_sol);
    for (int i = 1; i < data1_num; i++)
    {
        fprintf(gpid, "'%s' using 2:3 every :::%d::%d w l lw %f lc rgb 'black', " , data1.c_str(), i, i, lw_sol);
    }
    fprintf(gpid, "[%f:%f] x+%f w l lw %f lc rgb 'red' dt (16,8),",x_min, x_max, Pr_u-Pf_u, lw_vel); //　iso-velocity line
    fprintf(gpid, "'%s' using 1:2 every ::1::1 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=1
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=0
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 3.2 lc rgb 'white'\n",equi.c_str()); // equilibrium point h(t)=0

    // h(t)=0
    fprintf(gpid, "reset\n");
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/Phase0.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set size ratio '%f'\n", (y_max-y_min)/(x_max-x_min));
    fprintf(gpid, "set xrange [%f:%f]\n", x_min, x_max);
    fprintf(gpid, "set yrange [%f:%f]\n", y_min, y_max);
    fprintf(gpid, "set tics font '%s'\n", font_tics);
    fprintf(gpid,"set border back lw 2\n");
    fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 0.2 noborder behind\n", x_min, y_min, x_max, y_max);
    fprintf(gpid, "plot '%s' using 2:3 every :::%d::%d w l lw %f lc rgb 'black', " , data0.c_str(), 0, 0, lw_sol);
    for (int i = 1; i < data0_num; i++)
    {
        fprintf(gpid, "'%s' using 2:3 every :::%d::%d w l lw %f lc rgb 'black', " , data0.c_str(), i, i, lw_sol);
    }
    fprintf(gpid, "[%f:%f] x+%f w l lw %f lc rgb 'red' dt (16,8),",x_min, x_max, Pr_u-Pf_u, lw_vel); //　iso-velocity line
    fprintf(gpid, "'%s' using 1:2 every ::1::1 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=1
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=0
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 3.2 lc rgb 'white'\n",equi.c_str()); // equilibrium point h(t)=0

    fflush(gpid);
    pclose(gpid);
}