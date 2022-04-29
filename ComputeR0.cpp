#include <iostream>
#include <fstream>
//#include <filesystem>
#include <stdio.h> /* input, output, puts, NULL */
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <chrono>
#include <sstream>

using namespace std;


struct Vertex
{
    int Status_old = 0;
    int Status = 0;
    int Infector = -1;
    // Status: {S = 0 , I = 1, R=2}
    vector<int> adjList;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unifreal_dis(0.0, 1.0);

void CreateNetworkFromEdgeList(string FilePath, Vertex Nodes[]);
void CreateErdosReinyGraph(double p_grph, int NNodes, Vertex Nodes[]);
void SIRonNet(double beta, double mu, int NNodes,vector <int>& ListofNode, Vertex Nodes[]);
void ReadParameters(string FilePath,double parameters[]);
double mean(const std::vector<double> &v);
double sd(const std::vector<double> &v);

int main(int argc, char** argv)
{
    char label[] = "SIR";
    int NNode = 10000;

    Vertex Nodes[NNode];
    std::vector<int> ListofNode;
    for (int i=0; i<NNode; ++i) ListofNode.push_back(i); 

    double beta_s = 0.04 , beta_e = 0.05, beta_step = 0.1;
    double mu_s = 1.0/7.0, mu_e = (1.0/7.0)+0.01, mu_step = 0.2;
    int itr = 20;
    int NNodetoComputeR0 = 1000;

    double parameters[8] = { beta_s, beta_e, beta_step,
                            mu_s, mu_e, mu_step, 
                            (double)NNodetoComputeR0,(double)itr};

    string NetworkLabel = "Net";
    if (argc > 1) {
        string ConfigFilePath = argv[1];
        string NetworkFilePath = argv[2];
        CreateNetworkFromEdgeList(NetworkFilePath, Nodes);
        // Read parameters from file 
        ReadParameters(ConfigFilePath, parameters);
        beta_s = parameters[0], beta_e = parameters[1], beta_step = parameters[2] ;
        mu_s = parameters[3], mu_e = parameters[4], mu_step = parameters[5];
        NNodetoComputeR0 = (int)parameters[6];
        itr = (int)parameters[7];
        if (argc > 3){
            NetworkLabel =  argv[3];      
        }
    }
    
    double MeanDegree = 0; 
    for (int i = 0; i < NNode; i++)
    {
        MeanDegree += Nodes[i].adjList.size();
    }
    MeanDegree = MeanDegree / (double)NNode;

    
    cout << "parameters: "<<endl;
    cout << "beta= "<<beta_s<<":"<<beta_e<<":"<<beta_step<< endl;
    cout << "mu= "<<mu_s<<":"<<mu_e<<":"<<mu_step<< endl;
    cout << "Number Of nodes to Compute R0= "<<NNodetoComputeR0<<",  itr= " <<itr<<endl;

    //==============Header of the Result files ===============
    string HeaderFile = "beta,mu,node,R0_mean,R0_sd";

    //================Creat Result Files ====================
    //string Path = std::filesystem::current_path();
    char PathC[256];
    getcwd(PathC, 256);
    string Path = (string)PathC;
    char FileNameSuffix[128];
    snprintf(FileNameSuffix, sizeof(FileNameSuffix), "_k=%.4f.csv", MeanDegree );
    string FileName = "R0forNodes_" + NetworkLabel + (string)FileNameSuffix;
    string Filepath = Path + "/" + (string)FileName;
    ofstream file1;
    file1.open(Filepath);
    file1 << HeaderFile;
    file1 << endl;
    //========================================================
    //===========Creat Containers to store resuls=============
    std::vector<double> NofInfectedBynode_ItrVec(itr, 0);
    std::vector<std::array<double, 2>> R0nodeMeanSd_table;
    //========================================================
   
    vector<double> betaVec;
    vector<double> muVec;


    for (double betai = beta_s; betai < beta_e; betai += beta_step){betaVec.push_back(betai);}
    for (double mui = mu_s; mui <= mu_e; mui += mu_step){muVec.push_back(mui);}
  
    for (double beta : betaVec)
    {       
        auto start = chrono::steady_clock::now();
        for (double mu : muVec)
        {
            std::vector<std::array<double, 2>> R0nodeMeanSd_table;
            std::vector<int> ListofNode_sh = ListofNode;
            std::vector<int> ListofNodetoComputeR0;

            std::shuffle(ListofNode_sh.begin(), ListofNode_sh.end(), gen);
            int nii = 0;
            while (ListofNodetoComputeR0.size() < NNodetoComputeR0)
            {
                int nodei = ListofNode_sh[nii];
                if ((Nodes[nodei].adjList.size())>0){
                    ListofNodetoComputeR0.push_back(nodei);
                }
                nii+=1;
                // cout<<nii <<' '<< Nodes[nodei].adjList.size()<<endl;
            }
            
            for (int ni = 0; ni < NNodetoComputeR0; ni++)
            {
                int node = ListofNodetoComputeR0[ni];
                //auto start = chrono::steady_clock::now();
                for (int itrC = 0; itrC < itr; itrC++)
                {   
                    for (int i = 0; i < NNode; i++)
                    {
                        // Status: {S = 0, I = 1, R=2}
                        Nodes[i].Status = 0;
                        Nodes[i].Status_old = 0;
                        Nodes[i].Infector = -1;
                    }
                    Nodes[node].Status = 1;
                    Nodes[node].Status_old = 1;
                    int timestep = 0;
                    int NofInfectedBynode = 0;
                    while (Nodes[node].Status == 1 && timestep < 300)
                    {   
                        // Compute and recorde the current Status:
                    
                        // Update the old status with the current status and prepare to compute the next step.
                        for (size_t i = 0; i < NNode; i++)
                        {
                            Nodes[i].Status_old = Nodes[i].Status;  
                        }
                        timestep += 1;
                        // Compute the next status and update the current status:
                        SIRonNet(beta, mu, NNode, ListofNode, Nodes);
                    }
                    
                    for (int v : Nodes[node].adjList)
                    {
                        if(Nodes[v].Infector == node)
                        {
                            NofInfectedBynode +=1;
                        }
                    }
                    NofInfectedBynode_ItrVec[itrC] = NofInfectedBynode;
                    
                    
                }
                
                R0nodeMeanSd_table.push_back({mean(NofInfectedBynode_ItrVec),sd(NofInfectedBynode_ItrVec)});
                NofInfectedBynode_ItrVec.assign(NofInfectedBynode_ItrVec.size(),0);
          

            }
        
            for (int itt = 0; itt < R0nodeMeanSd_table.size(); ++itt)
            {
                file1 << beta << "," << mu<<",";
                file1 << ListofNodetoComputeR0[itt] ;
                for (int ittt = 0; ittt < R0nodeMeanSd_table[0].size(); ittt++)
                {
                    file1 << "," << R0nodeMeanSd_table[itt][ittt];
                }
                file1 << "\n";
            } 
            R0nodeMeanSd_table.clear();
            ListofNodetoComputeR0.clear();
            
        }
    
        
        auto end = chrono::steady_clock::now();
        auto diff = end - start;
        cout << "beta: "<< beta << ", run time: " << chrono::duration<double>(diff).count() << "s" << endl;

    }
          
    file1.close();
    cout << "The results stored at: " << Filepath << endl;
    
    return 0;
}

void CreateErdosReinyGraph(double p_grph, int NNodes, Vertex Nodes[])
{
    //std::uniform_real_distribution<double> unifreal_dis(0.0,1.0);
    for (int node1 = 0; node1 < NNodes; node1++)
    {

        for (int node2 = 0; node2 < node1; node2++)
        {
            double r = unifreal_dis(gen);
            if (r < p_grph)
            {
                Nodes[node1].adjList.push_back(node2);
                Nodes[node2].adjList.push_back(node1);
            }
        }
    }
}

void SIRonNet(double beta, double mu, int NNodes,vector <int>& ListofNode, Vertex Nodes[])
{
    // beta_k = P(S-->I) for kth starin
    // mu_k = P(I-->R) for kth starin
    for (int i : ListofNode)
    {
            if (Nodes[i].Status_old == 1) //if the source is in Status == I for kth strain (can be an Infector):
            {
                //Transmition
                std::vector<int> ListofNeighbour = Nodes[i].adjList;
                std::shuffle(ListofNeighbour.begin(), ListofNeighbour.end(), gen);
                for (int v : ListofNeighbour)
                {
                    if (Nodes[v].Status_old == 0) //if the target is in Status == S for kth strain:
                    {
                        double p_trans = beta;
                        double r1 = unifreal_dis(gen);
                        if (r1 < p_trans)
                        {
                            Nodes[v].Status = 1; //Status = I for kth strain
                            Nodes[v].Infector = i;
                        }
                    }
                }
                //Recovery
                double p_rec = mu;
                double r2 = unifreal_dis(gen);
                if (r2 < p_rec)
                {
                    Nodes[i].Status = 2; //Status = R for kth strain
                }
            }
    }
}



double mean(const std::vector<double> &v)
{
    double sum = 0;
    double len = v.size();
    for (auto &each: v)
        sum += each;

    return sum / len;
}

double sd(const std::vector<double> &v)
{
    double square_sum_of_difference = 0;
    double mean_var = mean(v);
    auto len = v.size();

    double tmp;
    for (auto &each: v) {
        tmp = each - mean_var;
        square_sum_of_difference += tmp * tmp;
    }

    return std::sqrt(square_sum_of_difference / (len - 1));
}

 void CreateNetworkFromEdgeList(string FilePath, Vertex Nodes[])
 {
    std::ifstream infile (FilePath);    // Load the file stream
    std::string line;                  // A line of values from text
    std::stringstream splitter;        // Prepare a stringstream as a splitter (splits on spaces) for reading key/values from a line
    
    // Make sure we can read the stream
    if (infile) {
        // As long as there are lines of data, we read the file
        while (std::getline(infile, line)) {
            int source, target;
            std::stringstream splitter;
            splitter << line;           // Load line into splitter
            //cout << line << endl ;           
            splitter >> source;         // Read the key back into temporary
            splitter >> target;         // Read the value back into temporary
            splitter.clear();           // Clear for next line
            // Add the edge to the Graph:
            Nodes[source].adjList.push_back(target);
            Nodes[target].adjList.push_back(source);
        
        }
    }
    else {
        // The file was not found or locked, etc...
        std::cout << "Unable to open file: " << FilePath << std::endl;
        exit(1);
    }
 }
    
void ReadParameters(string FilePath,double parameters[])
{
    std::ifstream infile (FilePath);    // Load the file stream
    std::string line;                  // A line of values from text
    std::stringstream splitter;        // Prepare a stringstream as a splitter (splits on spaces) for reading key/values from a line
    /*
    { beta_s, beta_e, beta_step,
                            mu_s, mu_e, mu_step, 
                            (double)NNodetoComputeR0,(double)itr};
    */
    // Make sure we can read the stream
    if (infile) {
        // As long as there are lines of data, we read the file
        while (std::getline(infile, line)) {
            string VariableName ;
            double tempDoublValue;
            splitter << line;           // Load line into splitter
            //cout << line << endl ;           
            splitter >> VariableName;         // Read the key back into temporary
           
            if(VariableName=="beta"){
                splitter >> tempDoublValue; 
                double beta_s = tempDoublValue ;
                splitter >> tempDoublValue; 
                double beta_e = tempDoublValue ;
                splitter >> tempDoublValue; 
                double beta_step = tempDoublValue ;
                parameters[0] = beta_s;
                parameters[1] = beta_e;
                parameters[2] = beta_step;
            } else if(VariableName=="mu"){
                splitter >> tempDoublValue; 
                double mu_s = tempDoublValue ;
                splitter >> tempDoublValue; 
                double mu_e = tempDoublValue ;
                splitter >> tempDoublValue; 
                double mu_step = tempDoublValue ;
                parameters[3] = mu_s;
                parameters[4] = mu_e;
                parameters[5] = mu_step;
            }else if (VariableName=="NNodetoComputeR0"){
                splitter >> tempDoublValue;   
                double NNodetoComputeR0 = tempDoublValue ;
                parameters[6] = NNodetoComputeR0;
            }else if (VariableName=="itr"){
                splitter >> tempDoublValue; 
                double itr = tempDoublValue ;
                parameters[7] = itr;
            }else {
                cout<<"Can not read all parameters from the config file. Please enter the parameters with the following keys:" << endl;
                cout<<"beta1, mu1, r2, tau2, r3, tau3, sigma3, t2, t3, Iinit, itr \n" << endl;
                break;
            }
            splitter.clear();           // Clear for next line   
        }

    }
    else {
        // The file was not found or locked, etc...
        std::cout << "Unable to open file: " << FilePath << std::endl;
        exit(1);
    }
}