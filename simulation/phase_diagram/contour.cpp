#include<string>
#include<stdlib.h>

using namespace std;

# define gr (double)1.0 

double vl_sec(double v, double gam);
double ac_sec(double a, double a_de, double a_de_s, double Dif, double rv, double gam, double Pf_s, double Pr_s);

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
    
    double Pf_s=para2[2], Pf_u=para2[0], Pr_s=para2[3], Pr_u=para2[1];
    double gf=para[0], gf_s=para[1], a_de=para[2], a_de_s=para[3], Dif=para[4], gam=para[5], rv=para[6];
    double x_min = Pf_u-0.1, x_max = Pf_s+0.1;
    double y_min = Pr_s-(Pf_s-Pf_u)+0.1, y_max = Pr_s+0.3;

    FILE *gpid;
    gpid = popen("gnuplot -persist","w"); // open gnuplot
    if(gpid == NULL)
    {
    printf("I can't find gnuplot.\n"); exit(1);
    }

    string equi = "result/param/param_plot.dat";

    // velocity
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/contour_vel.png'\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set border lw 2\n");
    fprintf(gpid, "set tics front out\n");
    fprintf(gpid, "set size ratio '%f'\n", (y_max-y_min)/(x_max-x_min));
    fprintf(gpid, "set xrange [%lf:%lf]\n", x_min, x_max);
    fprintf(gpid, "set yrange [%lf:%lf]\n", y_min, y_max);
    fprintf(gpid, "set xtics nomirror scale 1.5\n");
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 1 noborder behind\n", x_min, y_min, x_max, y_max);

    double vel, y_sec, lw_vel=5.;
    //　plot iso-velocity lines
    vel = -0.01;
    y_sec = vl_sec(vel, gam);
    fprintf(gpid, "plot [%f:%f] x+%f w l lw %lf lc rgb 'red' dt (5,2.5),",x_min, x_max, y_sec, lw_vel);
    vel = 0.0;
    y_sec = vl_sec(vel, gam);
    fprintf(gpid, "[%f:%f] x+%f w l lw %lf lc rgb 'red' dt (5,2.5),",x_min, x_max, y_sec, lw_vel);
    vel = 0.01;
    y_sec = vl_sec(vel, gam);
    fprintf(gpid, "[%f:%f] x+%f w l lw %lf lc rgb 'red' dt (5,2.5),",x_min, x_max, y_sec, lw_vel);
    vel = 0.02;
    y_sec = vl_sec(vel, gam);
    fprintf(gpid, "[%f:%f] x+%f w l lw %lf lc rgb 'red' dt (5,2.5),",x_min, x_max, y_sec, lw_vel);
    fprintf(gpid, "'%s' using 1:2 every ::1::1 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=１
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 4.2 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=0
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 3.0 lc rgb 'white'\n",equi.c_str()); // equilibrium point h(t)=0

    // acceleration
    fprintf(gpid, "set terminal pngcairo\n");
    fprintf(gpid, "set out 'result/fig/contour_acc.png'\n");
    fprintf(gpid, "reset\n");
    fprintf(gpid, "unset key\n");
    fprintf(gpid, "set border lw 2\n");
    fprintf(gpid, "set tics front out\n");
    fprintf(gpid, "set size ratio '%f'\n", (y_max-y_min)/(x_max-x_min));
    fprintf(gpid, "set xrange [%lf:%lf]\n", x_min, x_max);
    fprintf(gpid, "set yrange [%lf:%lf]\n", y_min, y_max);
    fprintf(gpid, "set xtics nomirror scale 1.5\n");
    fprintf(gpid, "set ytics nomirror scale 1.5\n");
    fprintf(gpid, "set object rect from %lf,%lf to %lf,%lf fillcolor rgb 'white' fillstyle solid 1 noborder behind\n", x_min, y_min, x_max, y_max); // solid：透明度

    double ac, slope=(a_de_s+Dif+rv*Dif)/(a_de+Dif+rv*Dif), lw_acc=5.;

    // plot iso-acceleration lines
    ac = -0.01;
    y_sec = ac_sec(ac, a_de, a_de_s, Dif, rv, gam, Pf_s, Pr_s);;
    fprintf(gpid, "plot [%f:%f] %lf*x+%f w l lw %lf lc rgb 'blue' dt (1.,1.),",x_min, x_max, slope, y_sec, lw_acc);
    ac = -0.005;
    y_sec = ac_sec(ac, a_de, a_de_s, Dif, rv, gam, Pf_s, Pr_s);;
    fprintf(gpid, "[%f:%f] %lf*x+%f w l lw %lf lc rgb 'blue' dt (1.,1.),",x_min, x_max, slope, y_sec, lw_acc);
    ac = 0.;
    y_sec = ac_sec(ac, a_de, a_de_s, Dif, rv, gam, Pf_s, Pr_s);;
    fprintf(gpid, "[%f:%f] %lf*x+%f w l lw %lf lc rgb 'blue' dt (1.,1.),",x_min, x_max, slope, y_sec, lw_acc);
    ac = 0.005;
    y_sec = ac_sec(ac, a_de, a_de_s, Dif, rv, gam, Pf_s, Pr_s);;
    fprintf(gpid, "[%f:%f] %lf*x+%f w l lw %lf lc rgb 'blue' dt (1.,1.),",x_min, x_max, slope, y_sec, lw_acc);
    fprintf(gpid, "'%s' using 1:2 every ::1::1 w p pt 7 ps 4.0 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=1
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 4.2 lc rgb 'black',",equi.c_str());  // equilibrium point h(t)=0
    fprintf(gpid, "'%s' using 1:2 every ::0::0 w p pt 7 ps 3.0 lc rgb 'white'\n",equi.c_str()); // equilibrium point h(t)=0
}

double vl_sec(double v, double gam)
// y section of iso-velocity line with velocity v
{
    return v/gam;
}

double ac_sec(double a, double a_de, double a_de_s, double Dif, double rv, double gam, double Pf_s, double Pr_s)
// y section of iso-acceleration line with acceleration a
{
    double c1=a_de+Dif+rv*Dif;
    double c2=-(a_de_s+Dif+rv*Dif)*Pf_s+(a_de+Dif+rv*Dif)*Pr_s;
    return (gam*c2-a)/(gam*c1);
}
