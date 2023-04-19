#include "msu_erhosampler/erhosampler.h"

using namespace NMSU_ERrhoSampler;

void CcorrVsEtaScott::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
	double Z=1.0/sqrt(2.0*PI);
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			corr(a,b)=0.0;
			if((a==0 ||a>3) && (b==0||b>3)){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==b){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==3 && (b==0||b>3)){
				corr(a,b)=Z*eta*exp(-eta*eta/2.0);
			}
			else if(b==3 && (a==0||a>3)){
				corr(a,b)=-Z*eta*exp(-eta*eta/2.0);
			}
		}
	}
}

void CcorrVsEtaOleh::init1(std::string path)
{
    std::cout<<"Begin reading input files:"<<std::endl;
    for(int i1=0;i1<this->Ncharges;i1++){
        std::cout<<this->charg_names[i1]<<" ";
        for(int j1=0;j1<this->Ncharges;j1++)
        {
            std::cout<<this->charg_names[j1]<<std::endl;
            std::string t_filename=path+charg_names[i1]+"-"+charg_names[j1]+".dat";
            std::cout<<t_filename;
            this->t_file.open(t_filename);
            std::cout<<" is open: "<< t_file.is_open()<<std::endl;
            for(int k=0;k<Netas;k++)
            {
                double t_1;
                t_file>>t_1;
                this->corr_input[Ncharges*i1+j1][k]=t_1;
            }
            this->t_file.close();
        }
    }
    
    
}
int CcorrVsEtaOleh::index_eta(double eta)
{
    double eta_max=10;
    double d_eta = 2*eta_max/(Netas-1);
    if (fabs(eta)>eta_max/2)
    {
        return Netas/4;
    }
    return (eta+eta_max)/d_eta;
}
int CcorrVsEtaOleh::index_charge(int a,int b)
{
    return Ncharges*a+b;
}
void CcorrVsEtaOleh::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
    
    int eta_number=index_eta(eta);
    for(a=0;a<7;a++){
        for(b=0;b<7;b++){
            corr(a,b)=0.0;
            if((a==0 ||a>3) && (b==0||b>3)){
                corr(a,b)=corr_input[index_charge(a,b)][eta_number];
            }
            else if(a==b){
                corr(a,b)=corr_input[index_charge(a,b)][eta_number];
            }
            else if(a==3 && (b==0||b>3)){
                corr(a,b)=corr_input[index_charge(a,b)][eta_number];
            }
            else if(b==3 && (a==0||a>3)){
                corr(a,b)=corr_input[index_charge(a,b)][eta_number];
            }
        }
        
    }

}
