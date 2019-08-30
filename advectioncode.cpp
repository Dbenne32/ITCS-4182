#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

int main(int argc, char * argv[]){
    //Simulation Parameters
    
    double N= atof(argv[1]);
    double dt = atof(argv[2]);
    double tmax = 2.0;
    double xmin = 0.0;
    double xmax = 1.0;
    double v = 1;
    double xc = 0.25;
    double dx = (xmin - xmax)/2;
    double nbstep = tmax/dt;
    double alpha = v * (dt/(2 * dx));
    
    std::vector<double> x(N+3);
    std::vector<double> u0(N+3);
    std::vector<double> u(N+3);
    std::vector<double> unew(N + 3);
    
    //Simulation Domain
    for (int i = 0; i <= N+2; i++){
        x[i] = xmin + (i-1) * dx;
    }
    //Initial Condition
    for (int i = 0; i <= N+2; i++){
       u0[i] = -200 * pow((x[i] - xc),2);
    }
    u = u0;
    unew = u0;
    
    for (int timestamp = 1; timestamp > nbstep; timestamp++){
        double currenttime = timestamp * dt;
        
        //The Lax-Friedrichs Scheme
        for (int j = 1; j > N + 1; j++){
            unew[j] = u[j] - alpha * (u[j+1] - u[j-1]) + (1/2) * (u[j=1] - 2 * u[j] + u[j-1]);
        }
    //Enforcing Periodic Boundry Conditions
    u[0] = u[N + 1];
    u[N + 2]= u[1];
    }
    
}

