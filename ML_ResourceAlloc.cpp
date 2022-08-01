#include <iostream>
#include <math.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include "twister.c"
using namespace std;

/*define Network Structure*/
#define M  10                    //the number of edges
#define N  5                     //the number of fogs for each eade
#define distance_FE 1000         //the distance between fog and edge (m)
#define distance_EC 100000       //the distance between edge and cloud (m)
#define training_B1 600
#define training_B2 100

/*define link queue id*/
#define FE 0        //link fog to edge
#define EC 1        //link edge to fog
#define CE 2        //link edge to cloud
#define EF 3        //link cloud to edge
#define FF 4        //link fog to fog
#define EE 5        //link edge to edge
#define F  6        //fog
#define E  7        //edge
#define C  8        //cloud
#define N_node 9

#define Predict 0
#define Train 1

/*define capacity*/
// server computation capacity (CPU cycles/s)
double upperbound_cloud = 100e+9;   //100GHz
double upperbound_edge = 15e+9;     //15GHz
double upperbound_fog = 2.5e+9;     //2.5GHz
// link computation capacity(bits/s)
double uplink_capacity = 1000000000;      //1Gps
double downlink_capacity = 1000000000;    //1Gps

double C_rsrc_T[N_node]={0.0};     // resource allocate for each node
double C_rsrc_P[N_node]={0.0};
double UB_node[N_node];             // set upperbound for each node
double LB_node_T[N_node];
double LB_node_P[N_node];

/*define system parameter*/
double lambda = 2000;                           //arrival rate of raw data
double mu_L1 = 15000.0;                        //average bits of raw data (bits)
double mu_L2 = 10560.0;                        //average bits of parameter  for model (bits)
double mu_T1 = 15000.0*training_B1*5;          //ML model training workload (CPU cycles)
double mu_P1 = 15000.0*1*0.5;                   //ML model inferencing workload (CPU cycles)
double mu_F1 = 10560.0*training_B2*0.1;          //Decentralized learning workload (CPU cycles)
double mu_F2 = 10560.0*training_B2*0.1;          //Federated Learning workload (CPU cycles)

/*define delay constraint*/
double constraint_T = 10.0;           //delay constraint of training (s)
double constraint_P = 0.05;         //delay constraint of inferencing (s)

/*define B1 B2*/
double B1 = 600.0;
double B2 = 100.0;

/*define workload arrival rate*/
double workload_P[N_node];
double workload_T[N_node][2];
double arrival_P[N_node];
double arrival_T[N_node][2];

/*define */
int scenario;
int numFed;

/*save to excel*/
double capacity[19][5];
double T_capacity[19][1];
double base_T[19][21] = {0.0};
double base_P[19][21] = {0.0};
double C_node_T[19][9] = {0.0};
double C_node_P[19][9] = {0.0};

// maximum capacity for each node
void init_max_rec(double M_fog, double M_edge, double M_cloud, double UL, double DL){
    UB_node[F] = M_fog;
    UB_node[E] = M_edge;
    UB_node[C] = M_cloud;
    UB_node[FE] = UL;
    UB_node[EC] = UL*N;
    UB_node[CE] = DL*N;
    UB_node[EF] = DL;
    UB_node[EE] = UL;
    UB_node[FF] = UL;
}
// set workload for each node
void set_workload(){

    for(int i=0; i<N_node; i++) workload_P[i] = 0.0;
    for(int i=0; i<N_node; i++){
        for(int j=0; j<2; j++){
            workload_T[i][j] = 0.0;
        }
    }
    int tag;
    tag = scenario % 3;

    //set prediction workload
    switch(tag){
        case 1:
            workload_P[FE] = mu_L1;
            workload_P[EC] = mu_L1;
            workload_P[C] = mu_P1;
            break;
        case 2:
            workload_P[FE] = mu_L1;
            workload_P[E] = mu_P1;
            break;
        case 0:
            workload_P[F] = mu_P1;
            break;
    }
    //set training workload
    switch(scenario){
        case 1: // NC/C
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L1;
            workload_T[C][0] = mu_T1;
            break;
        case 2: // NC/E
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L1;
            workload_T[CE][0] = mu_L2;
            workload_T[C][0] = mu_T1;
            break;
        case 3: // NC/F
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L1;
            workload_T[CE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[C][0] = mu_T1;
            break;
        case 4: // NE-C/C
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[E][0] = mu_T1;
            workload_T[C][0] = mu_F2;
            break;
        case 5: // NE-C/E
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[E][0] = mu_T1;
            workload_T[C][0] = mu_F2;
            break;
        case 6: // NE-C/F
            workload_T[FE][0] = mu_L1;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[E][0] = mu_T1;
            workload_T[C][0] = mu_F2;
            break;
        case 7: // NF-E-C/C
        case 8: // NF-E-C/E
        case 9: // NF-E-C/F
            workload_T[FE][0] = mu_L2;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[F][0] = mu_T1;
            workload_T[E][0] = mu_F2;
            workload_T[C][0] = mu_F2;
            break;
        case 10: // NF-C/C
        case 11: // NF-C/E
        case 12: // NF-C/F
            workload_T[FE][0] = mu_L2;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[F][0] = mu_T1;
            workload_T[C][0] = mu_F2;
            break;
        case 13:
            workload_T[FE][0] = mu_L2;
            workload_T[EC][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[EE][0] = mu_L2;
            workload_T[F][0] = mu_T1;
            workload_T[E][0] = mu_F1;
            workload_T[E][1] = mu_F2;
            break;
        case 14:
        case 15:
            workload_T[FE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[EE][0] = mu_L2;
            workload_T[F][0] = mu_T1;
            workload_T[E][0] = mu_F1;
            workload_T[E][1] = mu_F2;
            break;
        case 16:
        case 17:
        case 18:
            workload_T[FE][0] = mu_L2;
            workload_T[EC][0] = mu_L2;
            workload_T[CE][0] = mu_L2;
            workload_T[EF][0] = mu_L2;
            workload_T[FF][0] = mu_L2;
            workload_T[F][0] = mu_T1;
            workload_T[F][1] = mu_F1;
            workload_T[C][0] = mu_F2;

    }
}
// set arrival rate for each node
void set_arrival(){
    for(int i=0; i<N_node; i++) arrival_P[i] = 0.0;
    for(int i=0; i<N_node; i++){
        for(int j=0; j<2; j++){
            arrival_T[i][j] = 0.0;
        }
    }
    int tag;
    tag = scenario % 3;

    //set prediction arrival rate
    switch(tag){
        case 1:
            arrival_P[FE] = lambda/(M*N);
            arrival_P[EC] = lambda/M;
            arrival_P[C] = lambda;
            break;
        case 2:
            arrival_P[FE] = lambda/(M*N);
            arrival_P[E] = lambda/M;
            break;
        case 0:
            arrival_P[F] = lambda/(M*N);
            break;
    }
    //set training arrival rate
    switch(scenario){
        case 1: // NC/C
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/M;
            arrival_T[C][0] = lambda;
            break;
        case 2: // NC/E
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/M;
            arrival_T[CE][0] = lambda/B1;
            arrival_T[C][0] = lambda;
            break;
        case 3: // NC/F
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/M;
            arrival_T[CE][0] = lambda/B1;
            arrival_T[EF][0] = lambda/B1;
            arrival_T[C][0] = lambda;
            break;
        case 4: // NE-C/C
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/(M*B1);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[E][0] = lambda/M;
            arrival_T[C][0] = lambda/B1;
            break;
        case 5: // NE-C/E
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/(M*B1);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[E][0] = lambda/M;
            arrival_T[C][0] = lambda/B1;
            break;
        case 6: // NE-C/F
            arrival_T[FE][0] = lambda/(M*N);
            arrival_T[EC][0] = lambda/(M*B1);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[EF][0] = lambda/(B1*B2);
            arrival_T[E][0] = lambda/M;
            arrival_T[C][0] = lambda/B1;
            break;
        case 7: // NF-E-C/C
        case 8: // NF-E-C/E
        case 9: // NF-E-C/F
            arrival_T[FE][0] = lambda/(M*N*B1);
            arrival_T[EC][0] = lambda/(M*B1*B2);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[EF][0] = lambda/B1;
            arrival_T[F][0] = lambda/(M*N);
            arrival_T[E][0] = lambda/(M*B1);
            arrival_T[C][0] = lambda/(B1*B2);
            break;
        case 10: // NF-C/C
        case 11: // NF-C/E
        case 12: // NF-C/F
            arrival_T[FE][0] = lambda/(M*N*B1);
            arrival_T[EC][0] = lambda/(M*B1);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[EF][0] = lambda/(B1*B2);
            arrival_T[F][0] = lambda/(M*N);
            arrival_T[C][0] = lambda/B1;
            break;
        case 13:
            arrival_T[FE][0] = lambda/(M*N*B1);
            arrival_T[EC][0] = lambda/(B1*B2);
            arrival_T[EF][0] = lambda/(B1*B2);
            arrival_T[EE][0] = lambda/(M*B1*B2);
            arrival_T[F][0] = lambda/(M*N);
            arrival_T[E][0] = lambda/(M*B1);
            arrival_T[E][1] = lambda/(B1*B2);
            break;
        case 14:
        case 15:
            arrival_T[FE][0] = lambda/(M*N*B1);
            arrival_T[EF][0] = lambda/(B1*B2);
            arrival_T[EE][0] = lambda/(M*B1*B2);
            arrival_T[F][0] = lambda/(M*N);
            arrival_T[E][0] = lambda/(M*B1);
            arrival_T[E][1] = lambda/(B1*B2);
            break;
        case 16:
        case 17:
        case 18:
            arrival_T[FE][0] = lambda/(M*B2*B1);
            arrival_T[EC][0] = lambda/(M*B1*B2);
            arrival_T[CE][0] = lambda/(B1*B2);
            arrival_T[EF][0] = lambda/(B1*B2);
            arrival_T[FF][0] = lambda/(M*N*B1);
            arrival_T[F][0] = lambda/(M*N);
            arrival_T[F][1] = lambda/(M*B1);
            arrival_T[C][0] = lambda/(B1*B2);
            break;
    }
}
void set_lowerbound(){
    double alpha1;
    double alpha2;
    double unit = 1.0;
    //lowerbound for prediction
    for(int i=0; i<N_node; i++){
        LB_node_P[i] = arrival_P[i]*workload_P[i] + unit;
    }
    //set lowerbound for training
    for(int i=0; i<N_node; i++){
        if(workload_T[i][1] < 1e-15 && arrival_T[i][1] < 1e-15){
            LB_node_T[i] = arrival_T[i][0]*workload_T[i][0] + unit;
        }
        else{
            alpha1 = (arrival_T[i][0])/(arrival_T[i][0]+arrival_T[i][1]);
            alpha2 = (arrival_T[i][1])/(arrival_T[i][0]+arrival_T[i][1]);
            LB_node_T[i] = ((alpha1*workload_T[i][0])+(alpha2*workload_T[i][1]))*(arrival_T[i][0]+arrival_T[i][1]) + unit;
        }
    }
}
// calculate queue delay follow M/M/1
double queue_delay(double capacity, double demand, double arrival_rate){
    //printf("queue_delay:%lf\n",1/((capacity/demand) - arrival_rate));
    if( (1/((capacity/demand) - arrival_rate)) < 1e-15){
        printf("delay is less than zero\n");
        //exit(1);
    }
    return (1/((capacity/demand) - arrival_rate));
}
// calculate queue delay follow M/h2/1
double M_h2_1_delay(double arr1, double arr2, double dem1, double dem2, double capacity){
    double output = 0.0;
    double alpha1, alpha2;
    double mu_1, mu_2;
    double expected;
    double scquare_E;
    double utilization;

    //set alpha
    alpha1 = arr1/(arr1 + arr2);
    alpha2 = arr2/(arr1 + arr2);
    //set mu
    mu_1 = capacity/dem1;
    mu_2 = capacity/dem2;
    //calculate expected value 
    expected = (alpha1/mu_1) + (alpha2/mu_2);
    scquare_E = 2 * ((alpha1/(mu_1*mu_1)) + (alpha2/mu_2*mu_2));
    //calculate utilization (E[X]*lambda)
    utilization = expected * (arr1 + arr2);
    //printf("debug:%lf",utilization);
    //calculate delay
    output = expected + ((utilization*scquare_E)/(2*(1-utilization)*expected));
    //printf("M/H2/1:%lf\n",output);
    if(output < 1e-15){
        printf("delay is less than zero\n");
        //exit(1);
    }
    return output;
}
double propagation_delay(double distance){
    return (distance/((2.0/3.0)*3.0*100000000));
}
double waiting_time(double arrival_rate, int tag){

    //first-tier waiting time
    if(tag == 1){
        return (((B1-1)/2)*(1/arrival_rate));
    }
    //second-tier waiting time
    else if(tag == 2){
        return (((B1*B2)-1)/2*(1/arrival_rate));
    }
    else return 0;
}

// calculate inferencing delay
double inferencing_delay(double capacity[]){
    double communication_delay = 0.0;
    double computation_delay = 0.0;

    //communication_delay
    for(int i=FE; i<=EE; i++){
        if(workload_P[i] > 1e-15){
            communication_delay += queue_delay(capacity[i], workload_P[i], arrival_P[i]);
            if(i%2 == 0 && (i<=EF)) communication_delay += propagation_delay(distance_FE);
            else if(i%2 != 0  && (i<=EF)) communication_delay += propagation_delay(distance_EC);
        }
    }
    //printf("communication_delay: %lf\n",communication_delay);

    //computation delay
    for(int i=F; i<=C ;i++){
        if(workload_P[i] > 1e-15){
            computation_delay += queue_delay(capacity[i], workload_P[i], arrival_P[i]);
        }
    }
    //printf("computation_delay: %lf\n",computation_delay);
    return (communication_delay + computation_delay);

}

// calculate training delay
double training_delay(double capacity[]){
    double communication_delay = 0.0;
    double computation_delay = 0.0;

    //communication delay
    for(int i=FE; i<=EE; i++){
        if(workload_T[i][1] < 1e-15 && arrival_T[i][1] < 1e-15){ //check only one arrival in a queue
            if(workload_T[i][0] > 1e-15){
                communication_delay += queue_delay(capacity[i],workload_T[i][0], arrival_T[i][0]);
                if(i%2 == 0 && (i<=EF)) communication_delay += propagation_delay(distance_FE);
                else if(i%2 != 0  && (i<=EF)) communication_delay += propagation_delay(distance_EC);
            }
        } 
    }
    //computation delay
    for(int i=F; i<=C; i++){
        if(workload_T[i][1] < 1e-15 && arrival_T[i][1] < 1e-15){
            if(workload_T[i][0] > 1e-15){
                computation_delay += queue_delay(capacity[i],workload_T[i][0], arrival_T[i][0]);
            }
        }
        else{
            if(workload_T[i][0] > 1e-15){
                computation_delay += M_h2_1_delay(arrival_T[i][0], arrival_T[i][1], workload_T[i][0], workload_T[i][1], capacity[i]);
            }
        }
    }
    /*
    //if do first and second federation, need W1 and W2 
    if(numFed == 2){
        computation_delay += waiting_time(arrival_T[C][0],2);
        //printf("waiting_time: %lf\n",waiting_time(arrival_T[E][0],1)+waiting_time(arrival_T[C][0],2));
    }
    //do first federation, only need W1
    else if(numFed == 1){
        computation_delay += waiting_time(arrival_T[C][0],1);
        //printf("waiting_time: %lf\n",waiting_time(arrival_T[C][0],1));
    }*/
    return (communication_delay + computation_delay);
}

//average allocation(baseline method)
void average_alloc(double CR[], double MAX[], double wide, int todo){
    double UB[N_node], LB[N_node], interval[N_node], D_constraint, delay = 0.0;
    double arr[N_node] = {0.0};
    int flag;
    int passed_queue = 0;
    D_constraint = todo==Predict? constraint_P : constraint_T;
    for(int i=0; i<N_node; i++){
        UB[i] = MAX[i];
        if(todo == Predict){
            LB[i] = LB_node_P[i];
            arr[i] = arrival_P[i];
        }
        else if(todo == Train){
            LB[i] = LB_node_T[i];
            arr[i] = arrival_T[i][0];
        }
        if(arr[i]>1e-15){
            interval[i] = MAX[i] - LB[i];
            CR[i] = interval[i]*0.5+LB[i];
            passed_queue++;
        }
    }
    printf("passed_queue:%d\n",passed_queue);
    D_constraint = D_constraint/passed_queue;
    for(int i=0; i<N_node; i++){
        flag = 1;
        if(arr[i]>1e-15){
            if(todo == Predict){
                delay = queue_delay(CR[i], workload_P[i], arrival_P[i]);
                if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
            }
            else if(todo == Train){
                if(workload_T[i][1]>1e-15 && i>=F){
                    delay = M_h2_1_delay(arrival_T[i][0], arrival_T[i][1], workload_T[i][0], workload_T[i][1], CR[i]);
                }
                else{
                    delay = queue_delay(CR[i], workload_T[i][0], arrival_T[i][0]);
                    if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                    else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
                }
            }
            printf("delay_constraint: %.8f\n",D_constraint);
            while(delay > D_constraint || flag){
                //printf("1\n");
                if(delay <= D_constraint){
                    UB[i] = CR[i];
                    CR[i] = (CR[i] + LB[i])/2;
                }
                else{
                    LB[i] = CR[i];
                    CR[i] = (CR[i] + UB[i])/2;
                }
                if(todo == Predict){
                    delay = queue_delay(CR[i], workload_P[i], arrival_P[i]);
                    if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                    else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
                }
                else if(todo == Train){
                    if(workload_T[i][1]>1e-15 && i>=F){
                        delay = M_h2_1_delay(arrival_T[i][0], arrival_T[i][1], workload_T[i][0], workload_T[i][1], CR[i]);
                    }
                    else{
                        delay = queue_delay(CR[i], workload_T[i][0], arrival_T[i][0]);
                        if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                        else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
                    }
                }
                flag = 0;
                //printf("UB-LB:%.15f\n",((UB[i]-LB[i])/interval[i]));
                //printf("UB-LB-wide:%.15f\n",((UB[i]-LB[i])/interval[i])-wide);
                if(((UB[i]-LB[i])/interval[i]) - wide > 1e-15){
                    //printf("UB-LB:%.15f\n",(UB[i]-LB[i])/interval[i]);
                    //printf("wide:%.15f\n",wide);
                    //printf("%.15f\n",(UB[i]-LB[i])/interval[i] - wide);
                    flag = 1;
                }
                //printf("flag:%d\n",flag);
                //printf("delay:%8f\n",delay);
            }
        }
        else{
            continue;
        }

    }
    printf("\n");
    for(int i=0; i<N_node; i++){
        //printf("CR:%.8f\n",CR[i]);
        if(arr[i]>1e-15){
            if(todo == Predict){
                delay = queue_delay(CR[i], workload_P[i], arrival_P[i]);
                if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
                printf("resource:%.8f queue_delay:%.8f\n",CR[i],delay);
                base_P[scenario][i]=CR[i];
                base_P[scenario][i+N_node]=delay;
            }
            else if(todo == Train){
                if(workload_T[i][1]>1e-15 && i>=F){
                    delay = M_h2_1_delay(arrival_T[i][0], arrival_T[i][1], workload_T[i][0], workload_T[i][1], CR[i]);
                    printf("resource:%.8f queue_delay:%.8f\n",CR[i],delay);
                    base_T[scenario][i]=CR[i];
                    base_T[scenario][i+N_node]=delay;
                }
                else{
                    delay = queue_delay(CR[i], workload_T[i][0], arrival_T[i][0]);
                    if((i%2 == 0) && (i<=EF)) delay += propagation_delay(distance_FE);
                    else if((i%2 != 0) && (i<=EF)) delay += propagation_delay(distance_EC);
                    printf("resource:%.8f queue_delay:%.8f\n",CR[i],delay);
                    base_T[scenario][i]=CR[i];
                    base_T[scenario][i+N_node]=delay;
                }
            }
        }
    }
}

// cost function
double Cost(double current, double R, double MAX){
    double base=100.0;
    double cost_value;
    cost_value = pow(base, ((current+R)/MAX)-1) - pow(base, ((current)/MAX)-1);
    return cost_value;
}

// optimize bisection method
void optimize_bisection(double CR[], double MAX[], double wide, int todo){
    double UB[N_node], LB[N_node], interval[N_node], cost[N_node], D_constraint, delay;
    double arr[N_node] ={0.0};
    int flag=1;
    D_constraint = todo==Predict? constraint_P : constraint_T;
    for(int i=0; i<N_node; i++){
        UB[i] = MAX[i];
        if(todo == Predict){
            LB[i] = LB_node_P[i];
            arr[i] = arrival_P[i];
        }
        else if(todo == Train){
            LB[i] = LB_node_T[i];
            arr[i] = arrival_T[i][0];
            //printf("%.2f \n",LB[i]);
        }
        //initial: allocate resource between upperbound and lowerbound for each node
        if(arr[i]>1e-15){
            interval[i] = MAX[i] - LB[i];
            CR[i] = (MAX[i] + LB[i])/2;
        }
    }
    //printf("\n");
    //calculate delay
    delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
    int time=0;
    while(delay > D_constraint || flag){
        //calculate cost for each node
        for(int j=0; j<N_node; j++){
            if(arr[j]>1e-15){
                cost[j] = delay>D_constraint? Cost(CR[j], ((CR[j]+UB[j])/2)-CR[j], MAX[j]):Cost(CR[j], ((CR[j]+LB[j])/2)-CR[j], MAX[j]);
                //printf("capacity: %d %.8f\n",j,CR[j]);
            }
        }
        //choose the minimum cost of node
        int min_node;
        double min_cost;
        min_node = 999;
        min_cost = 999.0;
        for(int j=0; j<N_node; j++){
            if(arr[j]>1e-15){
                //printf("cost: %.16f\n",cost[j]);
                if(( (UB[j]-LB[j]) > wide) && (cost[j] < min_cost) ) {
                    min_node = j;
                    min_cost = cost[j];
                }
            }
        }
        //printf("min_cost:%.16f\n",min_cost);
        int adjust_node = min_node; //choose the minimum cost of node that passed
        //printf("adjust_node:%d\n",adjust_node);
        if(adjust_node == 999){
            double r[N_node]={0.0};
            double new_cost[N_node]={0.0};
            int mini_node = 999;
            double mini_cost = 999.0;
            for(int j=0; j<N_node; j++){
               if(arr[j]>1e-15){
                   r[j] = todo==Predict? workload_P[j]*0.5 : workload_T[j][0]*0.5;
                   //r[j] = 100;
                   new_cost[j] = Cost(CR[j],r[j], MAX[j]);
                   if(new_cost[j] < mini_cost){
                       mini_node = j;
                       mini_cost = new_cost[j];
                   }
               }
            }
            int choose_node = mini_node;
            //printf("adjust_node:%d\n",choose_node);
            CR[choose_node] += r[choose_node];
            delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
            //printf("delay:%.8f\n",delay);
            //printf("delay constraint is not satisfied");
            if(delay < D_constraint) break;
            else continue;
        }
        // choose the minimum cost node to do bisection
        if(delay <= D_constraint){
            UB[adjust_node] = CR[adjust_node];
            CR[adjust_node] = (CR[adjust_node] + LB[adjust_node])/2;
        }
        else{
            LB[adjust_node] = CR[adjust_node];
            CR[adjust_node] = (CR[adjust_node] + UB[adjust_node])/2;
        }
        //calculate the delay of updated capacity
        delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
        //printf("delay:%.8f\n",delay);
        //check the wide of upper bound and lower bound is close to bisection error for each node
        flag = 0;
        for(int i=0;i<N_node;i++){
            if(arr[i] > 1e-15){
                if((UB[i]-LB[i]) - wide > 1e-15){
                    flag = 1;
                }
            }
        }
        time ++;
    }
    printf("%d\n",time);
}

// optimize allocation method
void optimize_alloc(double CR[], double MAX[], int todo){
    double UB[N_node], LB[N_node], D_constraint, delay;
    double q_delay[N_node]={0.0};
    double cost[N_node] ={0.0};
    double arr[N_node] ={0.0};
    double R_unit = 100.0;
    //initialize
    D_constraint = todo==Predict? constraint_P : constraint_T;
    for(int i=0; i<N_node; i++){
        UB[i] = MAX[i];
        if(todo == Predict){
            LB[i] = LB_node_P[i];
            arr[i] = arrival_P[i];
        }
        else if(todo == Train){
            LB[i] = LB_node_T[i];
            arr[i] = arrival_T[i][0];
            //printf("%.2f \n",LB[i]);
        }
        //initial: allocate resource between upperbound and lowerbound for each node
        if(arr[i]>1e-15){
            CR[i] = LB[i];
        }
    }
    delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
    int time = 0;
    while(delay >= D_constraint){
        //calculate delay and cost of each node
        for(int j=0; j<N_node; j++){
            if(arr[j]>1e-15){
                //calculate delay of each node
                if(todo == Predict){
                    q_delay[j] = queue_delay(CR[j], workload_P[j], arrival_P[j])-queue_delay(CR[j]+R_unit, workload_P[j], arrival_P[j]);
                }
                else if(todo == Train){
                    if(workload_T[j][1]>1e-15 && j>=F){
                        q_delay[j] = M_h2_1_delay(arrival_T[j][0], arrival_T[j][1], workload_T[j][0], workload_T[j][1], CR[j])-M_h2_1_delay(arrival_T[j][0], arrival_T[j][1], workload_T[j][0], workload_T[j][1], CR[j]+R_unit);
                    }
                    else{
                        q_delay[j] = queue_delay(CR[j], workload_T[j][0], arrival_T[j][0])-queue_delay(CR[j]+R_unit, workload_T[j][0], arrival_T[j][0]);
                    }
                }
                //calculate cost of each node
                cost[j] = Cost(CR[j],R_unit,MAX[j]);
            }
        }
        //choose the maximum cost of node
        int max_node;
        double max_cost;
        max_node = -1;
        max_cost = -1.0;
        for(int j=0; j<N_node; j++){
            if(arr[j]>1e-15){
                //printf("delay/cost: %.16f\n",q_delay[j]/cost[j]);
                if(q_delay[j]/cost[j] > max_cost) {
                    max_node = j;
                    max_cost = q_delay[j]/cost[j];
                }
            }
        }
        //printf("max_cost:%.16f\n",max_cost);
        int adjust_node = max_node; //choose the maximum cost of node that passed
        //printf("adjust_node:%d\n",adjust_node);
        //update the capacity of adjusted node
        CR[adjust_node] += R_unit;
        //calclate the total delay
        delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
        //printf("delay:%.8f\n",delay);
        time++;
    }
    printf("times: %d\n",time);
}
void differ_bisection(double CR[], double MAX[], double wide, int todo){
    double UB[N_node], LB[N_node], interval[N_node], cost[N_node], D_constraint, delay;
    int flag=1;
    D_constraint = todo==Predict? constraint_P : constraint_T;
    for(int i=0; i<N_node; i++){
        UB[i] = MAX[i];
        if(todo == Predict){
            LB[i] = LB_node_P[i];
        }
        else if(todo == Train){
            LB[i] = LB_node_T[i];
            //printf("%.2f ",LB[i]);
        }
        interval[i] = MAX[i] - LB[i];
        CR[i] = interval[i]*genrand64()+LB[i];
        cost[i] = todo == Predict? 10*log(CR[i]-LB_node_P[i]) : 10*log(CR[i]-LB_node_T[i]);
    }
    //printf("\n");
    delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
    int time=0;
    while(delay > D_constraint || flag){
        //printf("delay:%.8f\n",delay);
        int min_node,max_node,min_cost,max_cost;
        min_node = max_node = 999;
        min_cost = 999;
        max_cost = 0;
        for(int j=0; j<N_node; j++){
            //printf("cost:%.8f \n",cost[j]);
            if(( (((UB[j]-LB[j])/interval[j]) > wide) && cost[j] < min_cost) ){
                min_node = j;
                min_cost = cost[j];
            }
            if((((UB[j]-LB[j])/interval[j]) > wide) && cost[j] > max_cost){
                max_node = j;
                max_cost = cost[j];
            }
        }
        int adjust_node = (delay <= D_constraint)? max_node : min_node;
        //printf("adjust_node:%d\n",adjust_node);
        if(adjust_node == 999) adjust_node = rand()%(C-0+1)+0;
        if(delay <= D_constraint){
            UB[adjust_node] = CR[adjust_node];
            CR[adjust_node] = (CR[adjust_node] + LB[adjust_node])/2;
        }
        else{
            LB[adjust_node] = CR[adjust_node];
            CR[adjust_node] = (CR[adjust_node] + UB[adjust_node])/2;
        }
        cost[adjust_node] = todo == Predict? 10*log(CR[adjust_node]-LB_node_P[adjust_node]) : 10*log(CR[adjust_node]-LB_node_T[adjust_node]);
        delay = todo==Predict? inferencing_delay(CR) : training_delay(CR);
        flag = 0;
        for(int i=0;i<N_node;i++) if(((UB[i]-LB[i])/interval[i]) - wide > 1e-15) flag = 1;
        time ++;
    }
    printf("%d\n",time);
}
void report(double CR_P[], double CR_T[]){
    double T_delay=0.0, P_delay=0.0;
    int training_node;
    double totoal_weighted_capacity=0.0;
    double W_FE=2.0, W_EF=2.0, W_EC=2.0, W_CE=2.0, W_FF=1.0, W_EE=1.0, W_C=4.0, W_E=8.0, W_F=12.0;
    /*
    for(int i=0; i<N_node; i++){
        if(workload_T[i][0] < 1e-15){
            CR_T[i] = 0.0;
        }
        if(workload_P[i] < 1e-15){
            CR_P[i] = 0.0;
        }
    }*/
    for(int i=0; i<N_node; i++){
        C_node_T[scenario][i] = CR_T[i];
        C_node_P[scenario][i] = CR_P[i];
    }
    P_delay = inferencing_delay(CR_P);
    T_delay = training_delay(CR_T);
    printf("=== Scenario %d ===\n",scenario);
    printf("All link capacity\n");
    for(int i=0; i<N_node; i++) printf("%.8f ",CR_T[i]);
    printf("\n");
    for(int i=0; i<N_node; i++) printf("%.8f ",CR_P[i]);
    printf("\n");
    printf("Capacity of Fog: %.8f %.8f\n",CR_P[F],CR_T[F]);
    printf("Capacity of Edge: %.8f %.8f\n",CR_P[E],CR_T[E]);
    printf("Capacity of Cloud: %.8f %.8f\n",CR_P[C],CR_T[C]);
    printf("Total Capacity: %.8f\n", (CR_P[C]+CR_T[C]) + M*(CR_P[E]+CR_T[E]) + M*N*(CR_P[F]+CR_T[F]));
    totoal_weighted_capacity = W_C*(CR_P[C]+CR_T[C]) + M*W_E*(CR_P[E]+CR_T[E]) + M*N*W_F*(CR_P[F]+CR_T[F]) + M*N*(W_FE*(CR_P[FE]+CR_T[FE])+W_EF*(CR_P[EF]+CR_T[EF])) + M*(W_EC*(CR_P[EC]+CR_T[EC])+W_CE*(CR_P[CE]+CR_T[CE])) + (M*(M-1)/2)*W_FF*(CR_P[FF]+CR_T[FF]) + (N*(N-1)/2)*W_EE*(CR_P[EE]+CR_T[EE]);
    printf("Total Weighted Capacity: %.8f\n",totoal_weighted_capacity);
    printf("Training Delay: %.16f\n",T_delay);
    printf("Prediction Delay: %.16f\n\n",P_delay);

    capacity[scenario][0] = CR_P[F] + CR_T[F];
    capacity[scenario][1] = CR_P[E] + CR_T[E];
    capacity[scenario][2] = CR_P[C] + CR_T[C];
    capacity[scenario][3] = (CR_P[C]+CR_T[C]) + M*(CR_P[E]+CR_T[E]) + M*N*(CR_P[F]+CR_T[F]);
    capacity[scenario][4] = totoal_weighted_capacity;
    if(scenario == 1 || scenario == 2 || scenario == 3){
        training_node = C;
        T_capacity[scenario][0] = CR_T[training_node];
    }
    else if(scenario == 4 || scenario == 5 || scenario == 6){
        training_node = E;
        T_capacity[scenario][0] = CR_T[training_node]*M;
    }
    else{
        training_node = F;
        T_capacity[scenario][0] = CR_T[training_node]*M*N;
    }
}
int main(){
    
    int choose;
    printf("=================================\n");
    printf("Choose Algorithm\n");
    printf("(1) Base Allocation\n");
    printf("(2) Cost Bisection\n");
    printf("(3) Optimize Allocation\n");
    printf("=================================\n");
    printf("Enter number: ");
    scanf("%d",&choose);
    int ID;
    double zero[N_node];
    double bisection_error = 100;
    double delay_error = 0.00001;
    int training_node;

    FILE *fp1;
    FILE *fp2;
    FILE *fp3;
    switch(choose){
        case 1:
            fp1 = fopen("base_capacity.csv","w");
            fp2 = fopen("average_C_node_P.csv","w");
            fp3 = fopen("average_C_node_T.csv","w");
            for(ID=1; ID<=18; ID++){
                init_genrand64(880703);
                scenario = ID;
                //initialize
                set_workload();
                set_arrival();
                init_max_rec(upperbound_fog, upperbound_edge, upperbound_cloud, uplink_capacity, downlink_capacity);
                set_lowerbound();
                for(int i=0; i<N_node; i++) C_rsrc_P[i] = C_rsrc_T[i] = zero[i] = 0.0;
                printf("Lower_bound\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_T[i]);
                }
                printf("\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_P[i]);
                }
                printf("\n");
                average_alloc(C_rsrc_P, UB_node, 0.0000001, Predict);
                average_alloc(C_rsrc_T, UB_node, 0.0000001, Train);
                report(C_rsrc_P,C_rsrc_T);
            }
            break;
        case 2:
            fp1 = fopen("bisection_capacity.csv","w");
            fp2 = fopen("bisection_C_node_P.csv","w");
            fp3 = fopen("bisection_C_node_T.csv","w");
            for(ID=1; ID<=18; ID++){
                init_genrand64(880703);
                scenario = ID;
                //initialize
                set_workload();
                set_arrival();
                init_max_rec(upperbound_fog, upperbound_edge, upperbound_cloud, uplink_capacity, downlink_capacity);
                set_lowerbound();
                for(int i=0; i<N_node; i++) C_rsrc_P[i] = C_rsrc_T[i] = zero[i] = 0.0;
                printf("Lower_bound\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_T[i]);
                }
                printf("\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_P[i]);
                }
                printf("\n");
                //bisection search algorithm
                optimize_bisection(C_rsrc_P, UB_node, bisection_error, Predict);
                optimize_bisection(C_rsrc_T, UB_node, bisection_error, Train);
                //differ_bisection(C_rsrc_P, UB_node, bisection_error, Predict)
                //differ_bisection(C_rsrc_T, UB_node, bisection_error, Train);
                report(C_rsrc_P,C_rsrc_T);
            }
            break;
        case 3:
            fp1 = fopen("optimizeAlloc_capacity.csv","w");
            fp2 = fopen("optimizeAlloc_C_node_p.csv","w");
            fp3 = fopen("optimizeAlloc_C_node_T.csv","w");
            for(ID=1; ID<=18; ID++){
                init_genrand64(880703);
                scenario = ID;
                //initialize
                set_workload();
                set_arrival();
                init_max_rec(upperbound_fog, upperbound_edge, upperbound_cloud, uplink_capacity, downlink_capacity);
                set_lowerbound();
                for(int i=0; i<N_node; i++) C_rsrc_P[i] = C_rsrc_T[i] = zero[i] = 0.0;
                printf("Lower_bound\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_T[i]);
                }
                printf("\n");
                for(int i=0; i<N_node; i++){
                    printf("%.8f ",LB_node_P[i]);
                }
                printf("\n");
                optimize_alloc(C_rsrc_P, UB_node, Predict);
                optimize_alloc(C_rsrc_T, UB_node, Train);
                report(C_rsrc_P,C_rsrc_T);
            }
            break;
        default:
            printf("Enter number: ");
            scanf("%d",&choose);
    }
    fprintf(fp1,"%s,%s,%s,%s,%s,%s\n","ID","Fog Capacity","Edge Capacity", "Cloud Capacity", "Total Capacity", "Total Weighted Capacity");
    fprintf(fp2,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","ID","Node FE","Node EC", "Node CE", "Node EF", "Node FF", "Node EE", "Node F", "Node E", "Node C");
    fprintf(fp3,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","ID","Node FE","Node EC", "Node CE", "Node EF", "Node FF", "Node EE", "Node F", "Node E", "Node C");
    for(int i=1; i<=18; i++){
        fprintf(fp1,"%d, %.8f,%.8f,%.8f,%.8f,%.8f\n",i,capacity[i][0],capacity[i][1],capacity[i][2],capacity[i][3],capacity[i][4]);
        fprintf(fp2,"%d, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n",i, C_node_P[i][0],C_node_P[i][1],C_node_P[i][2],C_node_P[i][3],C_node_P[i][4],C_node_P[i][5],C_node_P[i][6],C_node_P[i][7],C_node_P[i][8]);
        fprintf(fp3,"%d, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n",i, C_node_T[i][0],C_node_T[i][1],C_node_T[i][2],C_node_T[i][3],C_node_T[i][4],C_node_T[i][5],C_node_T[i][6],C_node_T[i][7],C_node_T[i][8]);
    }
}