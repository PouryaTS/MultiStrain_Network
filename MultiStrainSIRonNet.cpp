#include <iostream>
#include <fstream>
#include <filesystem>
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

using namespace std;

#define NStrain 3
#define NState (int)pow(3, NStrain)

struct Vertex
{
    int Status_old[NStrain] = {0};
    int Status[NStrain] = {0};
    int Infector = -1;
    // Status: {S = 0 , I = 1, R=2}
    vector<int> adjList;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> unifreal_dis(0.0, 1.0);

void CreateNetworkFromEdgeList(string FilePath, Vertex Nodes[]);
void CreateErdosReinyGraph(double p_grph, int NNodes, Vertex Nodes[]);
void InitializingSeeds(int NNodes, int Nstrains, int Nseeds[], Vertex Nodes[], const bool ResetNodes = true);
void InitializingSeeds2(int NNodes, int Nstrains, int Nseeds[], Vertex Nodes[], const bool ResetNodes = true);
void MultiStrainSIRonNet(double beta[], double mu[], double sigma[][NStrain], int NNodes,vector <int>& ListofNode, Vertex Nodes[]);
int MapState2DecimalNumber(int State[], int Nstrains);
int myPow(int x, int p);

int main(int argc, char** argv)
{
    char label[] = "SIR";
    int NNode = 10000;

    Vertex Nodes[NNode];
    std::vector<int> ListofNode;
    for (int i=0; i<NNode; ++i) ListofNode.push_back(i); // 1 2 3 4 5 6 7 8 9

    double MeanDeg = 5;
    double p_grph = (double)MeanDeg / (double)NNode; // p = 0.005
    CreateErdosReinyGraph(p_grph, NNode, Nodes);
    double R0_1 = 1.5, mu_1 = 1.0/7.0, tau2 = 1,tau3 = 1, r3 = 1;
    //double r2 = 1.8;
    double r2_s = 1.0 , r2_e = 2.0;
    if (argc > 1) {
    r2_s = std::atof(argv[1]);
    r2_e = std::atof(argv[2]);
    }

    int Nseeds[NStrain] = {50, 0, 0};

   
    const bool ProduceEventMatix = false;

    //==================================================================================================
    //  - Map the States to a decimal number Index.
    //  - Map the States Index to the Column Index of file.
    //  - Extract Valid States (Infection coexistance is forbbiden) and Infected States (SSI,SIS,...).
    //  - Create labels of the states and the header of files.
    //==================================================================================================

    std::array<string, NState> LabelArray;         //Label of each state (S->0, I->1, R->2)
    int MapState2Index[3][3][3];                   // Maping Matrix-Maps states to a decimal index (000=0, 001 = 1,002 = 2, 010=3, ...).
    std::vector<int> ValidStatus;                  //Valid States (Infection coexistance is forbbiden)
    std::vector<int> InfectedStatus;               //Infected States (SSI,SIS,...)
    std::array<int, NState> MapStateIndx2ColIndx1; // Index1 -> Column index of pupulation in each state
    std::array<int, NState> MapStateIndx2ColIndx2; // Index2 -> Column index of Infection Incident for each Infected state.

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                string Lab;
                Lab += label[i];
                Lab += label[j];
                Lab += label[k];
                int nI = 0;
                (i == 1) ? (nI += 1) : (nI += 0);
                (j == 1) ? (nI += 1) : (nI += 0);
                (k == 1) ? (nI += 1) : (nI += 0);
                int temp[3] = {i, j, k};
                int ind = MapState2DecimalNumber(temp, 3);
                MapState2Index[i][j][k] = ind;
                LabelArray[ind] = Lab;

                if (nI < 2) // Valid States (No coexistance of Infection)
                {
                    ValidStatus.push_back(ind);
                    if (nI == 1) // Infected Status
                    {
                        InfectedStatus.push_back(ind);
                    }
                }
            }
        }
    }

    MapStateIndx2ColIndx1.fill(-1);
    for (int it = 0; it < ValidStatus.size(); it++)
    {
        MapStateIndx2ColIndx1[ValidStatus[it]] = it + 2;
    }
    int maxIndx1 = *std::max_element(MapStateIndx2ColIndx1.begin(), MapStateIndx2ColIndx1.end());
    MapStateIndx2ColIndx2.fill(-1);
    for (int it = 0; it < InfectedStatus.size(); it++)
    {
        MapStateIndx2ColIndx2[InfectedStatus[it]] = it + maxIndx1 + 1;
    }

    //==============Header of the Result files ===============
    string HeaderFile1 = "itr,t,node,status_previous,status_current,Infector,status_Infector";
    string HeaderFile2;
    HeaderFile2 += "r2,Sigma,t1,t2,it,t";
    for (int it = 0; it < ValidStatus.size(); it++)
    {
        HeaderFile2 += ",";
        HeaderFile2 += ("N_" + LabelArray[ValidStatus[it]]);
    }
    for (int it = 0; it < InfectedStatus.size(); it++)
    {
        HeaderFile2 += ",";
        HeaderFile2 += ("Inc_" + LabelArray[InfectedStatus[it]]);
    }

    //================Creat Result Files ====================
    string Path = std::filesystem::current_path();
    char FileName1[128];
    snprintf(FileName1, sizeof(FileName1), "TransmitionTrack_3Strain_R1=%.2f_mu1=%.2f_tau=%.2f_r2=%.1f:%0.1f.csv", R0_1, mu_1, tau2, r2_s, r2_e);
    string Filepath1 = Path + "/" + (string)FileName1;
    ofstream file1;
    if (ProduceEventMatix)
    {
        file1.open(Filepath1);
        file1 << HeaderFile1;
        file1 << endl;
    }
    char FileName2[128];
    snprintf(FileName2, sizeof(FileName2), "TimeSerie_3Strain_R1=%.2f_mu1=%.2f_tau=%.2f_r2=%.1f:%0.1f.csv", R0_1, mu_1, tau2, r2_s, r2_e);
    string Filepath2 = Path + "/" + (string)FileName2;
    ofstream file2;
    file2.open(Filepath2);
    file2 << HeaderFile2;
    file2 << endl;
    //========================================================
    //===========Creat Containers to store resuls=============
    int maxIndx2 = *std::max_element(MapStateIndx2ColIndx2.begin(), MapStateIndx2ColIndx2.end());
    const int StateArraySize = maxIndx2;
    std::vector<int> State_current((StateArraySize + 1), 0);
    std::vector<std::vector<int>> Res_timeserie_table;
    std::array<int, 7> Res_TransmitionTrack;
    std::vector<std::array<int, 7>> Res_TransmitionTrack_table;
    //========================================================
    int IndxSSS = MapStateIndx2ColIndx1[MapState2Index[0][0][0]];
    int IndxISS = MapStateIndx2ColIndx1[MapState2Index[1][0][0]];
    int IndxSIS = MapStateIndx2ColIndx1[MapState2Index[0][1][0]];
    int IndxSSI = MapStateIndx2ColIndx1[MapState2Index[0][0][1]];
    int IndxRSS = MapStateIndx2ColIndx1[MapState2Index[2][0][0]];

    
    vector<double> rVec;
    vector<double> sigmaVec;
    for (double ri = r2_s; ri < r2_e; ri += 0.1){rVec.push_back(ri);}
    for (double sigmai = 0.0; sigmai < 1.001; sigmai += 0.1){sigmaVec.push_back(sigmai);}
    
    for (double r2 : rVec)
    {       
        for (double Sigma3 : sigmaVec)
        {
    double mu_2 = mu_1 / tau2;
    double mu_3 = mu_1 / tau3;
    double beta_1 = R0_1 * mu_1 / MeanDeg;
    double beta_2 = (r2 * R0_1) * mu_2 / MeanDeg;
    double beta_3 = (r3 * R0_1) * mu_3 / MeanDeg;

    //double R0_12 = MeanDeg * (beta_1 + beta_2) / (mu_1 + mu_2);
    double beta[NStrain] = {beta_1, beta_2, beta_3};
    double mu[NStrain] = {mu_1, mu_2, mu_3};
    double Sigma12 = 0.05;
    //ouble Sigma3 = 0.5;
    double Sigma[3][NStrain] = {
        {Sigma12, Sigma12, Sigma3},
        {Sigma12, Sigma12, Sigma3},
        {Sigma3, Sigma3, Sigma3}};
    
    auto start = chrono::steady_clock::now();
    int itr = 100;
    for (int t1 = 0; t1 < 150; t1 += 15)
    {
        for (int t2 = t1; t2 < (t1+150); t2 += 15)
        {
            //auto start = chrono::steady_clock::now();
            for (int itrC = 0; itrC < itr; itrC++)
            {
                InitializingSeeds2(NNode, NStrain, Nseeds, Nodes, true);
                int timestep = 0;
                int NofInfc = 0;
                for (int k = 0; k < NStrain; k++)
                {
                    NofInfc = NofInfc + Nseeds[k];
                }

                int flag_emerge2 = 0;
                int flag_emerge3 = 0;
                int SSS0 = NNode - NofInfc;
                while (NofInfc > 0 && timestep < 1000)
                {   
                    // Emerge new strains
                    int SSS1 = State_current[IndxSSS];
                    if (timestep == t1 && flag_emerge2 == 0)                   
                    {
                        SSS0 = SSS1;
                        flag_emerge2 = 1;
                        int NseedNewVariant[NStrain] = {0, 50, 0};
                        InitializingSeeds2(NNode, NStrain, NseedNewVariant, Nodes, false);
                        //cout << timestep << " flag2: " << R0_1 << " " << SSS0 << endl;
                    }
                    //if (((double)SSS1 / (double)SSS0) < (1 / (pt2 * R0_12)) && flag_emerge2 == 1 && flag_emerge3 == 0)
                    if (timestep == t2 && flag_emerge3 == 0)
                    {
                        flag_emerge3 = 1;
                        int NseedNewVariant[NStrain] = {0, 0, 50};
                        InitializingSeeds2(NNode, NStrain, NseedNewVariant, Nodes, false);
                        //cout << timestep<< " flag3: "<<R0_12<<" "<< SSS1 <<" "<<SSS0<< endl;
                    }   

                    // Compute and recorde the current Status:
                    State_current.assign(State_current.size(), 0);
                    State_current[0] = itrC;
                    State_current[1] = timestep;
                    for (int i = 0; i < NNode; i++)
                    {
                        int Status_o = MapState2Index[Nodes[i].Status_old[0]][Nodes[i].Status_old[1]][Nodes[i].Status_old[2]];
                        int Status_n = MapState2Index[Nodes[i].Status[0]][Nodes[i].Status[1]][Nodes[i].Status[2]];
                        bool StatusChange = (Status_o != Status_n);
                        double alpha_Ik = 1;
                        for (int kk = 0; kk < NStrain; kk++)
                        {
                            if (Nodes[i].Status[kk] == 1)
                            {
                                alpha_Ik = 0;
                                break;
                            }
                        }
                        bool existsI = (alpha_Ik == 0);
                        if (StatusChange && existsI)
                        {
                            int Infector = Nodes[i].Infector;
                            int Status_Infctor = MapState2Index[Nodes[Infector].Status_old[0]][Nodes[Infector].Status_old[1]][Nodes[Infector].Status_old[2]];
                            // Recording the transmition track
                            if (ProduceEventMatix)
                            {
                                Res_TransmitionTrack = {itrC, timestep, i, Status_o, Status_n, Infector, Status_Infctor};
                                Res_TransmitionTrack_table.push_back(Res_TransmitionTrack);
                            }
                            int Ind2 = MapStateIndx2ColIndx2[Status_n];
                            State_current[Ind2] += 1;
                        }
                        int Ind1 = MapStateIndx2ColIndx1[Status_n];
                        State_current[Ind1] += 1;
                    }
                    Res_timeserie_table.push_back(State_current);

                    // Update the old status with the current status and prepare to compute the next step.
                    for (size_t i = 0; i < NNode; i++)
                    {
                        for (int k = 0; k < NStrain; k++)
                        {
                            Nodes[i].Status_old[k] = Nodes[i].Status[k];
                        }
                    }

                    NofInfc = 0;
                    for (int s : InfectedStatus)
                    {
                        NofInfc += State_current[MapStateIndx2ColIndx1[s]];
                    }
                    timestep += 1;

                    // Compute the next status and update the current status:
                    MultiStrainSIRonNet(beta, mu, Sigma, NNode, ListofNode, Nodes);
                }
            }

            if (ProduceEventMatix)
            {
                for (int itt = 0; itt < Res_TransmitionTrack_table.size(); ++itt)
                {
                    file1 << r2 << "," << Sigma3<<",";
                    file1 << t1 << "," <<t2;
                    for (int ittt = 0; ittt < Res_TransmitionTrack_table[0].size(); ittt++)
                    {
                        file1 << "," << Res_TransmitionTrack_table[itt][ittt];
                    }
                    file1 << "\n";
                }
            }
            for (int itt = 0; itt < Res_timeserie_table.size(); ++itt)
            {
                file2 << r2 << "," << Sigma3<<",";
                file2 << t1 << "," <<t2;
                for (int ittt = 0; ittt < Res_timeserie_table[0].size(); ittt++)
                {
                    file2 << "," << Res_timeserie_table[itt][ittt];
                }
                file2 << "\n";
            }

            Res_timeserie_table.clear();
            Res_TransmitionTrack_table.clear();

            /*auto end = chrono::steady_clock::now();
            auto diff = end - start;
            cout << "t1: "<< t1 << ", t2: " << t2 << ", run time: " << chrono::duration<double>(diff).count() << "s" << endl;*/
        }
    }
            auto end = chrono::steady_clock::now();
            auto diff = end - start;
            cout << "r2: "<< r2 << ", Sigma: " << Sigma3 << ", run time: " << chrono::duration<double>(diff).count() << "s" << endl;

        }
    }


    cout << ProduceEventMatix << endl;
    if (ProduceEventMatix)
    {
        file1.close();
        cout << "The results stored at: " << Filepath1 << endl;
    }
    file2.close();
    cout << "The results stored at: " << Filepath2 << endl;

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

void MultiStrainSIRonNet(double beta[], double mu[], double sigma[][NStrain], int NNodes,vector <int>& ListofNode, Vertex Nodes[])
{
    // beta_k = P(S-->I) for kth starin
    // mu_k = P(I-->R) for kth starin
    std::random_shuffle(ListofNode.begin(), ListofNode.end());

    for (int i : ListofNode)
    {
        for (int k = 0; k < NStrain; k++)
        {
            if (Nodes[i].Status_old[k] == 1) //if the source is in Status == I for kth strain (can be an Infector):
            {
                //Transmition
                for (int v : Nodes[i].adjList)
                {
                    if (Nodes[v].Status_old[k] == 0) //if the target is in Status == S for kth strain:
                    {
                        double alpha_Ik = 1;
                        for (int kk = 0; kk < NStrain; kk++)
                        {
                            if (Nodes[v].Status[kk] == 1)
                            {
                                alpha_Ik = 0;
                                break;
                            }
                        }
                        int StateIndx_v = 9 * Nodes[v].Status_old[0] + 3 * Nodes[v].Status_old[1] + Nodes[v].Status_old[2];
                        double alpha_Rk = 0;
                        if (StateIndx_v == 0)
                        {
                            alpha_Rk = 1;
                        }
                        else if (StateIndx_v == 18 || StateIndx_v == 8)
                        {
                            alpha_Rk = sigma[0][k];
                        }
                        else if (StateIndx_v == 6 || StateIndx_v == 20)
                        {
                            alpha_Rk = sigma[1][k];
                        }
                        else if (StateIndx_v == 2 || StateIndx_v == 24)
                        {
                            alpha_Rk = sigma[2][k];
                        }

                        double p_trans = alpha_Ik * alpha_Rk * beta[k];
                        double r1 = unifreal_dis(gen);
                        if (r1 < p_trans)
                        {
                            Nodes[v].Status[k] = 1; //Status = I for kth strain
                            Nodes[v].Infector = i;
                        }
                    }
                }
                //Recovery
                double p_rec = mu[k];
                double r2 = unifreal_dis(gen);
                if (r2 < p_rec)
                {
                    Nodes[i].Status[k] = 2; //Status = R for kth strain
                }
            }
        }
    }
}

int MapState2DecimalNumber(int State[], int Nstrains)
{
    int sum = 0;
    for (int k = 0; k < Nstrains; k++)
    {
        sum += State[k] * myPow(3, (Nstrains - 1 - k));
    }
    return sum;
}
int myPow(int x, int p)
{
    int i = 1;
    for (int j = 1; j <= p; j++)
        i *= x;
    return i;
}

void InitializingSeeds(int NNodes, int Nstrains, int Nseeds[], Vertex Nodes[], const bool ResetNodes)
{

    std::uniform_int_distribution<int> unifint_dis(0, (NNodes-1));

    vector<int> vec;
    int NseedT = 0;
    for (int k = 0; k < Nstrains; k++)
    {
        NseedT = NseedT + Nseeds[k];
    }

    int n = 0;
    while (n < NseedT)
    {
        int number = unifint_dis(gen);
        auto result1 = std::find(begin(vec), end(vec), number);
        if (result1 == std::end(vec))
        {
            vec.push_back(number);
            n++;
        }
    }
    if (ResetNodes)
    {
        for (int i = 0; i < NNodes; i++)
        {
            for (int k = 0; k < Nstrains; k++)
            {
                Nodes[i].Status[k] = 0;
                Nodes[i].Status_old[k] = 0;
                // Status: {S = 0, I = 1, R=2}
            }
            Nodes[i].Infector = -1;
        }
    }

    int startpoint = 0;
    int endpoint = 0;
    for (int k = 0; k < Nstrains; k++)
    {
        endpoint = startpoint + Nseeds[k];
        if (startpoint != endpoint)
        {
            std::vector<int> split_i(vec.begin() + startpoint, vec.begin() + endpoint);
            for (int v : split_i)
            {
                Nodes[v].Status[k] = 1;
                Nodes[v].Status_old[k] = 1;
                // Status: {S = 0, S = 1, R=2}
            }
        }
        startpoint = endpoint;
    }
}

void InitializingSeeds2(int NNodes, int Nstrains, int Nseeds[], Vertex Nodes[], const bool ResetNodes)
{
    std::uniform_int_distribution<int> unifint_dis(0, (NNodes-1));

    if (ResetNodes)
    {
        for (int i = 0; i < NNodes; i++)
        {
            for (int k = 0; k < Nstrains; k++)
            {
                Nodes[i].Status[k] = 0;
                Nodes[i].Status_old[k] = 0;
                // Status: {S = 0, I = 1, R=2}
            }
            Nodes[i].Infector = -1;
        }
    }

    for (int k = 0; k < Nstrains; k++)
    {
        int Nseedi = Nseeds[k];
        int n = 0;
        while (n < Nseedi)
        {
            int number = unifint_dis(gen);
            double alpha_Ik = 1;
            for (int kk = 0; kk < NStrain; kk++)
            {
                if (Nodes[number].Status[kk] == 1 || Nodes[number].Status_old[kk] == 1)
                {
                    alpha_Ik = 0;
                    break;
                }
            }
            if (Nodes[number].Status[k] ==0 && alpha_Ik == 1)
            {
                Nodes[number].Status[k] = 1;
                Nodes[number].Status_old[k] = 1;
                // Status: {S = 0, S = 1, R=2}
                n++;
            }
        }
    }
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
    