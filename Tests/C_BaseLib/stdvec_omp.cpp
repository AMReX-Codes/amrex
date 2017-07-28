#include<iostream>
#include<vector>
int main()
{
    int n=10;
    std::vector<double> usrdat;
    std::vector<int> globalkey;

    usrdat.resize(n);

    for(int i=0;i<n;i++)
    {
        usrdat[i]=std::rand()%1000;
    }
    
    std::cout<<"usrdat:";
    for(int i=0;i<n;i++)
    {
        std::cout<<usrdat[i]<<"\t";
    }
    std::cout<<"\n";

    std::vector<double> globalvec;

#pragma omp parallel
    {
        std::vector<double> localvec;
        std::vector<double> localkey;

#pragma omp for
        for(int i=0;i<n;i++)
        {
            localvec.push_back(usrdat[i]);
            localkey.push_back(i);
        }

#pragma omp critical
        {
            globalvec.insert(globalvec.end(),localvec.begin(),localvec.end());
            globalkey.insert(globalkey.end(),localkey.begin(),localkey.end());
        }
    }

    std::cout<<"globalvec:";
    for(int i=0;i<n;i++)
    {
        std::cout<<globalvec[globalkey[i]]<<"\t";
    }
    std::cout<<"\n";
}
