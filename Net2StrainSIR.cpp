#include <iostream>
#include <fstream>
#include <stdio.h>      /* input, output, puts, NULL */
#include <math.h> 
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
//#include <EigenRand/EigenRand>
#include <Eigen/Dense>




using namespace std;

std::random_device rd;
std::mt19937 gen(rd());

struct Vertex
{
    int Status_old = 0; 
    int Status = 0;
    int Infector = -1;
    // Status: {S = 0 , I_s = 1, I_f=2, R_s=3, R_f=4, I_f_s=5, I_s_f=6, R=7}
    vector<int> adjList;
};

std::uniform_real_distribution<double> unifreal_dis(0.0,1.0);
void CreateErdosReinyGraph(double p_grph,int NNodes,Vertex Nodes[]);
void InitializingNodes(int NNodes,int Nseed_1,int Nseed_2, Vertex Nodes[]);
void Net2StrainSIR(double beta_f, double beta_s, double mu_f, double mu_s, double sigma, int NNodes, Vertex Nodes[]);
//template <size_t rows, size_t cols>
//void Writ2DArr2csv(std::string filename, double (&array)[rows][cols]);
void Writ2DArr2csv2(std::string filename, double** array,int rows, int cols);


int main(){

    
    int NNode = 10000;
    Vertex Nodes[NNode];
    int Nseed1 = 50, Nseed2 = 50;

    double MeanDeg = 5;
    double p_grph = (double)MeanDeg/(double)NNode; //0.005
    CreateErdosReinyGraph(p_grph, NNode, Nodes);
    
    double R0_f = 1.8, mu_f = 0.6, tau = 1.5;
    vector<double> rVec;
    vector<double> sigmaVec;
    for (double ri = 0.6; ri <= 1.6; ri+=0.1){rVec.push_back(ri);}
    for (double sigmai = 0.0; sigmai <= 1; sigmai+=0.1){sigmaVec.push_back(sigmai);}
    int rows = rVec.size();
    int cols = sigmaVec.size();
    
    double **ResMatrix1;
    ResMatrix1 = new double *[rows];
    for(int i = 0; i <cols; i++)
        ResMatrix1[i] = new double[cols];

    double **ResMatrix2;
    ResMatrix2 = new double *[rows];
    for(int i = 0; i <cols; i++)
        ResMatrix2[i] = new double[cols];

    int rC=0;
    int sigmaC = 0;
    
    string Path = string(get_current_dir_name());
    string fileTransTrack = Path+"/TransmitionTrack.csv";
    string fileTimeSerie = Path+"/TimeSerie.csv";
    fstream file1(fileTransTrack,ios::out);
    fstream file2(fileTimeSerie,ios::out);
for (double r:rVec)
{   
    
    int sigmaC = 0;
    for (double sigma:sigmaVec)
   {  

    
    double beta_f = R0_f*mu_f;
    double mu_s = mu_f/tau;
    double beta_s = (r*R0_f)*mu_s;
    
    int itr = 2;
    int IcidenceMatrix[itr][2]={0}; 
    std::array<int, 10> Res_timeserie;
    std::vector<std::array<int, 10>> Res_timeserie_table;
    std::array<int, 7> Res_TransmitionTrack;
    std::vector<std::array<int, 7>> Res_TransmitionTrack_table;

    for (int itrC = 0; itrC < itr; itrC++){
    
        InitializingNodes(NNode, Nseed1, Nseed2, Nodes);
        int timestep = 0;
        int NofInfc = Nseed1+Nseed2;
        int State_old[8] = {NNode-NofInfc, Nseed1, Nseed2,0,0,0,0,0};
    
        int Incidence_s = Nseed1;
        int Incidence_f = Nseed2;
    
        while (NofInfc>0 && timestep<1000)
        {
            // Compute the next status and update the current status:
            Net2StrainSIR(beta_f, beta_s, mu_f, mu_s, sigma, NNode,Nodes);
            
            int State_new[8] = {0};
            for (int i = 0; i < NNode; i++) {
                if((Nodes[i].Status_old != Nodes[i].Status) && (Nodes[i].Status == 1||Nodes[i].Status == 2||Nodes[i].Status ==5||Nodes[i].Status ==6)){
                    int Status_o = Nodes[i].Status_old ;
                    int Status_n = Nodes[i].Status ;
                    int Infector = Nodes[i].Infector ;
                    int Status_Infctor = Nodes[Infector].Status_old ;
                    Res_TransmitionTrack = {itrC,timestep,i,Status_o ,Status_n, Infector,Status_Infctor};
                    Res_TransmitionTrack_table.push_back(Res_TransmitionTrack);
                }
                
                int S = Nodes[i].Status;
                State_new[S] +=1;
            }
            Res_timeserie[0]=itrC;
            Res_timeserie[1]=timestep;
            for (int i = 2; i < 10; i++){Res_timeserie[i]=State_new[i-2];}
            
            Res_timeserie_table.push_back(Res_timeserie);

            // Update the old status with the current status and prepare to compute the next step.
            for (size_t i = 0; i < NNode; i++) {Nodes[i].Status_old = Nodes[i].Status;}
            
            int New_s = ((State_new[1]-State_old[1])<0) ? 0 : (State_new[1]-State_old[1]);
            int New_fs = ((State_new[5]-State_old[5])<0) ? 0 : (State_new[5]-State_old[5]);
            int IncidenceNew_s = New_s+New_fs;

            int New_f = ((State_new[2]-State_old[2])<0) ? 0 : (State_new[2]-State_old[2]);
            int New_sf = ((State_new[6]-State_old[6])<0) ? 0 : (State_new[6]-State_old[6]);
            int IncidenceNew_f = New_f+New_sf;

            for (int i = 0; i < 8; i++){State_old[i]=State_new[i];}
            timestep +=1;
            NofInfc =State_new[1]+State_new[2]+State_new[5]+State_new[6];

            Incidence_s = Incidence_s + IncidenceNew_s;
            Incidence_f = Incidence_f + IncidenceNew_f;
    
        /*    cout << "t: "<<timestep << " state: ";
            for (int i = 0; i < 8; i++){cout<<State_new[i]<<" ";}
            cout <<"Inc_f:" <<Incidence_f << " Inc_s:" << Incidence_s<< endl;
        */

        }
        
        IcidenceMatrix[itrC][0] = Incidence_s;
        IcidenceMatrix[itrC][1] = Incidence_f;

        //cout << itrC<<" t:"<< timestep <<" Inc_f:" <<Incidence_f << " Inc_s:" << Incidence_s<< endl;
    }
    /*
    for (int i = 0; i < itr; i++){
        cout<<"itr: "<<i<<" Inc_s:"<< IcidenceMatrix[i][0]<<" Inc_f:"<<IcidenceMatrix[i][1]<<endl;
        }
    */
    double Sum1 =0,Sum2=0;
    for (int i = 0; i < itr; i++)
    {
        Sum1+=IcidenceMatrix[i][0];
        Sum2+=IcidenceMatrix[i][1];  
    }
    Sum1 = ((double)Sum1) / (double)itr;
    Sum2 = ((double)Sum2) / (double)itr;

    ResMatrix1[rC][sigmaC] = Sum1;
    ResMatrix2[rC][sigmaC] = Sum2;

    for ( int itt = 0 ; itt < Res_TransmitionTrack_table.size(); ++itt)
    {
        file1 << r <<","<< sigma;
        for (int ittt = 0; ittt < 7; ittt++)
        { 
            file1 << ","<<Res_TransmitionTrack_table[itt][ittt];
        }
        file1 << "\n";
    }
    for ( int itt = 0 ; itt < Res_timeserie_table.size(); ++itt)
    {
        file2 << r << "," << sigma;
        for (int ittt = 0; ittt < 10; ittt++)
        {
            file2 << ","<<Res_timeserie_table[itt][ittt];
        }
        file2 << "\n";
    }
/*
    for ( int itt = 0 ; itt < Res_timeserie_table.size(); ++itt)
    {
        for (int ittt = 0; ittt < 10; ittt++)
        {
        cout << Res_timeserie_table[itt][ittt]<< ",";
        }
        cout << endl;
    }*/
    

        sigmaC+=1;
        cout << "r="<<r<<" sigma="<<sigma<<endl;
   }
   rC+=1;
   
}

file1.close();
file2.close();  

    for(int i=0;i<rows;i++){
        cout << ResMatrix1[i][0];
        for(int j=1;j<cols;j++){
            cout<<", "<<ResMatrix1[i][j];   
        }  
        cout<<"\n"; 
    }

    //string Path = string(get_current_dir_name());
    char FileName1[128];
    snprintf(FileName1,sizeof(FileName1), "SIR2Strain_Inc_S_N=%d.csv",NNode);
    string Filepath1 = Path+"/"+(string)FileName1;
    char FileName2[128];
    snprintf(FileName2,sizeof(FileName2), "SIR2Strain_Inc_F_N=%d_Rf=%f_muf=%f_tau=%f.csv",NNode,R0_f,mu_f,tau);
    string Filepath2 = Path+"/"+(string)FileName2;

    Writ2DArr2csv2(Filepath1, ResMatrix1,rows, cols);
    Writ2DArr2csv2(Filepath2, ResMatrix2,rows, cols);
    
    for (int i = 0; i < rows; i++) {
        delete [] ResMatrix1[i];}

    delete [] ResMatrix1;
    
    for (int i = 0; i < rows; i++) {
        delete [] ResMatrix2[i];}

    delete [] ResMatrix2;

    return 0;   

}


void CreateErdosReinyGraph(double p_grph,int NNodes,Vertex Nodes[]){
    
    //std::uniform_real_distribution<double> unifreal_dis(0.0,1.0);
    for (int node1 = 0; node1 < NNodes; node1++){

        for (int node2 = 0; node2 < node1; node2++){
            double r = unifreal_dis(gen);
            if (r < p_grph)
            {
              Nodes[node1].adjList.push_back(node2);
              Nodes[node2].adjList.push_back(node1);   
            }
            
        }
        
    }
}

void InitializingNodes(int NNodes,int Nseed_1,int Nseed_2, Vertex Nodes[]){
 
    std::uniform_int_distribution<int> unifint_dis(0,NNodes);

    vector<int> vec; 
    int Nseed = Nseed_1+Nseed_2;
    int n=0;
    while (n<Nseed){
        int number = unifint_dis(gen);
        auto result1 = std::find(begin(vec), end(vec), number);
        if (result1 == std::end(vec))
        {
            vec.push_back(number);
            n++;
        } 
    }
    std::vector<int> split_1(vec.begin(), vec.begin() + Nseed_1);
    std::vector<int> split_2( vec.begin() + Nseed_1,vec.end());
    //cout<< split_lo[0]<<endl;
    //split_1[0]=0;
    //split_2[0]=1;
    for (int i = 0; i < NNodes; i++)
    {
        Nodes[i].Status = 0;
        Nodes[i].Status_old = 0;
        // Status: {S = 0, I_s = 1, I_f=2}
    }
    
    for (int v: split_1)
    {
        Nodes[v].Status = 1;
        Nodes[v].Status_old = 1;
        // Status: {S = 0, I_s = 1, I_f=2}
    }
    for (int v: split_2)
    {
        Nodes[v].Status = 2;    
        Nodes[v].Status_old = 2;
        // Status: {S = 0, I_s = 1, I_f=2}
    }
}

void Net2StrainSIR(double beta_f, double beta_s, double mu_f, double mu_s, double sigma, int NNodes, Vertex Nodes[]){
 
    int NofNodesWithInfc = 0;
    // beta_f = P(S-->If)
    // beta_s = P(S-->Is)
    double pt_s_f = sigma * beta_f ;       // P(Rs-->Isf)
    double pt_f_s = sigma * beta_s ;       // P(Rf-->Ifs)
    // mu_f = P(If-->Rf) and P(Isf-->R)
    // mu_s = P(Is-->Rs) and P(Ifs-->R)

    //Vertex NodeTemp[NNodes];
    //for (size_t i = 0; i < NNodes; i++) {NodeTemp[i] = Nodes[i];}


    for (size_t i = 0; i < NNodes; i++){
        
        if (Nodes[i].Status_old == 1 || Nodes[i].Status_old == 5 ){ // if Status == I_s or I_fs
            // transmition 
            for (int v: Nodes[i].adjList) {
                double r1 = unifreal_dis(gen);
                if (Nodes[v].Status_old == 0 && r1 < beta_s ){ // if Status==S
                    Nodes[v].Status = 1; //Status = I_s
                    Nodes[v].Infector = i;
                    }
                if (Nodes[v].Status_old == 4 && r1 < pt_f_s ){ // if Status==R_f
                    Nodes[v].Status = 5; // Status = I_fs
                    Nodes[v].Infector = i;
                    }
            }
            // recovery
            double r2 = unifreal_dis(gen);
            if(Nodes[i].Status_old == 1 && r2 < mu_s  ){ // if Status == I_s
                Nodes[i].Status = 3; // Status = R_s
            } else if(Nodes[i].Status_old == 5 && r2 < mu_s ){ // if Status == I_fs
                Nodes[i].Status = 7; // Status = R
            }
            
        } else if (Nodes[i].Status_old == 2 || Nodes[i].Status_old == 6 ){ // if Status == I_f or I_sf
            // transmition
            for (int v: Nodes[i].adjList) {
                double r1 = unifreal_dis(gen);
                if (Nodes[v].Status_old == 0 && r1 < beta_f  ){ // if Status==S
                    Nodes[v].Status = 2; //Status = I_f
                    Nodes[v].Infector = i;
                    }
                if (Nodes[v].Status_old == 3 && r1 < pt_s_f  ){ // if Status==R_s
                    Nodes[v].Status = 6; // Status = I_sf
                    Nodes[v].Infector = i;
                    }
            }
            // recovery
            double r2 = unifreal_dis(gen);
            if(Nodes[i].Status_old == 2 && r2 < mu_f ){// if Status == I_f
                Nodes[i].Status = 4; // Status = R_f
            } else if(Nodes[i].Status_old == 6 && r2 < mu_f ){ // if Status == I_sf
                Nodes[i].Status = 7; // Status = R
            }
        } 
    }
    
}

/*
template <size_t rows, size_t cols>
void Writ2DArr2csv(std::string filename, double (&array)[rows][cols]){
    fstream file(filename,ios::out);
    
    for(int i=0;i<rows;i++){
        file << array[i][0];
        for(int j=1;j<cols;j++){
            file<<", "<<array[i][j];   
        }  
        file<<"\n"; 
    }
    file.close();  
    cout<<"The results stored at: "<<filename<<endl;
}
*/

void Writ2DArr2csv2(std::string filename, double** array,int rows, int cols){
    fstream file(filename,ios::out);
    
    for(int i=0;i<rows;i++){
        file << array[i][0];
        for(int j=1;j<cols;j++){
            file<<", "<<array[i][j];   
        }  
        file<<"\n"; 
    }
    file.close();  
    cout<<"The results stored at: "<<filename<<endl;
}