//This program is to perform HiC compartmentalization in raw matrix file
//the process is to fit nij=Ci*Cj
//Only inter-chr reads are used. i.e. results from each chr are all connected, but no distance relation considered. Rabl could potentially have an effect.
//Incorporates arbitrary n compartments, and combine B and C back to one parameter, simplifying the calculation.
//Incorporates the effect of Rabl explicitly.
//Rabl structure is explicitly parameterized as R=1+aRiRj+b*(Ri^2Rj+RiRj^2), where Ri, Rj each take value between -0.5 (Centromere) and 0.5 (Telomere) .
//Incorporates translocated region exclusion. All interactions between regions in the same line in the exlcuded.txt are excluded from the calculation.
//added removal ENCODE blacklisted regions

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <math.h>
#include <omp.h>
#include "twister.h"
using namespace std;

const double maxc=0.9999;
const int thr_interaction=1;
static int nRabl=10;
const double steplength=0.04;
static int Nsub;

typedef vector<vector<double> > Matrix;

void parseStr(string& s, char c, vector<string>& vs)
{
	vs.clear();
	for (int i=0;;++i)
	{
		string::size_type pos=s.find(c);
		if (pos==string::npos)	
		{
			vs.push_back(s);
			break;
		}
		else
		{
			string s1=s.substr(0,pos);
			vs.push_back(s1);
			string s2=s.substr(pos+1,s.size()-pos-1);
			s=s2;
		}
	}
}

double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void LUdecomp(const Matrix& mat, Matrix& md)
{
	int n=mat.size();
	md=mat;
	for (int j=0;j<n;++j)
	{
		for (int i=0;i<=j;++i)
		{
			double s=0;
			for (int k=0;k<i;++k)
			{
				s+=md[i][k]*md[k][j];
			}
			md[i][j]=mat[i][j]-s;
		}
		for (int i=j+1;i<n;++i)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[i][k]*md[k][j];
			}
			md[i][j]=(mat[i][j]-s)/md[j][j];
		}
	}
}


//|m|
double det(const Matrix& mat)
{
	if (mat.size()==0)
	{
		return 1;
		cerr<<"Not a Matrix!"<<endl;
	}
	else if (mat.size()==1)
	{
		return mat[0][0];
	} 
	else if (mat.size()==2)
	{
		return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
	}
	else if (mat.size()==3)
	{
		return mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+mat[0][2]*mat[1][0]*mat[2][1]
			-(mat[0][2]*mat[1][1]*mat[2][0]+mat[0][1]*mat[1][0]*mat[2][2]+mat[0][0]*mat[1][2]*mat[2][1]);
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=1.0;
		for (int i=0;i<mat.size();++i)
		{
			s*=md[i][i];
		}
		return s;
	}
}

double innerprod(const vector<double>& x)
{
	int n=x.size();
	double d=0;
	for(int i=0;i<n;++i)
	{
		d+=x[i]*x[i];	
	}
	return d;
}

double innerprod(const vector<double>& x1, const vector<double>& x2)
{
	int n=x1.size();
	int n2=x2.size();
	if (n!=n2)
	{
		cerr<<"unmatched dimension in inner production!"<<endl;
		return 0;	
	}
	double d=0;
	for(int i=0;i<n;++i)
	{
		d+=x1[i]*x2[i];	
	}
	return d;
}

//m^-1
/*
Matrix inverse(const Matrix& mat)
{
	Matrix m1;
	m1=mat;
	if (mat.size()<=0) 
	{
		cerr<<"Not a Matrix!"<<endl;
	}
	Matrix md(mat);
	LUdecomp(mat,md);
	int n=mat.size();
	for (int j=0;j<n;++j)
	{
		double s=0;
		for (int k=0;k<j;++k)
		{
			s+=md[j][k]
		}
	}
}*/

//(m^-1)x

vector<double> invp(const Matrix& mat, const vector<double>& x)
{
	if (mat.size()!=x.size() || mat.size()==0)
	{
		cerr<<"Unmatched dimensions for scalar production!"<<endl;
		return x;
	}
	else if (mat.size()==1)
	{
		vector<double> d(1,x[0]/mat[0][0]);
		return d;
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=0;
		int n=mat.size();
		vector<double> c(x);
		for (int j=0;j<n;++j)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[j][k]*c[k];
			}
			c[j]=x[j]-s;
		}
		double sum=0;
		vector<double> d(x);
		for (int j=n-1;j>=0;--j)
		{
			double s=0;
			for (int k=j+1;k<n;++k)
			{
				s+=md[j][k]*d[k];
			}
			d[j]=(c[j]-s)/md[j][j];
		}
		return d;
	}
}

double scpro(const Matrix& mat, const vector<double>& x)
{
	if (mat.size()!=x.size() || mat.size()==0)
	{
		cerr<<"Unmatched dimensions for scalar production!"<<endl;
		return 0;
	}
	else if (mat.size()==1)
	{
		return x[0]*x[0]/mat[0][0];
	}
	else if (mat.size()==2)
	{
		return (mat[1][1]*x[0]*x[0]-(mat[1][0]+mat[0][1])*x[0]*x[1]+mat[0][0]*x[1]*x[1])/det(mat);
	}
	else if (mat.size()==3)
	{
		return ((mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])*x[0]*x[0]
			-(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0]+mat[0][1]*mat[2][2]-mat[0][2]*mat[2][1])*x[0]*x[1]
			+(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]+mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1])*x[0]*x[2]
			+(mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0])*x[1]*x[1]
			-(mat[0][0]*mat[1][2]-mat[0][2]*mat[1][0]+mat[0][0]*mat[2][1]-mat[0][1]*mat[2][0])*x[1]*x[2]
			+(mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])*x[2]*x[2])
			/det(mat);
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=0;
		int n=mat.size();
		vector<double> c(x);
		for (int j=0;j<n;++j)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[j][k]*c[k];
			}
			c[j]=x[j]-s;
		}
		double sum=0;
		vector<double> d(x);
		for (int j=n-1;j>=0;--j)
		{
			double s=0;
			for (int k=j+1;k<n;++k)
			{
				s+=md[j][k]*d[j];
			}
			d[j]=(c[j]-s)/md[j][j];
			sum+=x[j]*d[j];
		}
		return sum;
	}
}

void raninit(vector<double>& x)
{
	int n=x.size();
	vector<double> v(n);
	double s=0;
	for (int i=0;i<n;++i)
	{
		v[i]=-log(1-rndu());
		s+=v[i];
	}	
	for (int i=0;i<n;++i)
	{
	  x[i]=v[i]/s;
	}
}


//variable-length window
class VLW
{
public:
	string chrName;
	int chrStart;
	int chrEnd;
	int vlwsnumber;
	double Rabl;
	int Rabl_index;
	int color; //0:R,1:O,2:Y,3:G,4:B,5:P
	bool null;
	bool blackListed;
	int Excluded;
};

class Block
{
public:
	string chrName;
	int chrStart;
	int chrEnd;
	int blockStart; //[blockStart,blockEnd) is the region to be blocked
	int blockEnd;	
};

bool leftNonOverlap(const VLW& loc1, const VLW& loc2)
{
	if (loc1.chrName<loc2.chrName) return true;
	else if (loc1.chrName==loc2.chrName && loc1.chrEnd<=loc2.chrStart) return true;
	return false;
}

void readVLWs(const string& ifname, vector<VLW>& vlws)
{
	vlws.clear();
	ifstream fin(ifname.c_str());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<4) continue;
		VLW vlw1;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str());
		vlw1.chrEnd=atoi(vs[2].c_str());
		if (vs[3]=="RED") vlw1.color=0;
		if (vs[3]=="ORA") vlw1.color=1;
		if (vs[3]=="YEL") vlw1.color=2;
		if (vs[3]=="GRE") vlw1.color=3;
		if (vs[3]=="BLU") vlw1.color=4;
		if (vs[3]=="PUR") vlw1.color=5;
		vlws.push_back(vlw1);
	}
}

int readVLWindow(const string& ifname, vector<VLW>& vlws)
{
	vlws.clear();
	ifstream fin(ifname.c_str());
	if (!fin.is_open())
	{
		cout<<"Bad window File Name"<<endl;
		return -1;	
	}
	string s;
	while (!fin.eof())
	{
		getline(fin,s);
		if (s==""||s.substr(0,3)!="chr") continue;
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<4) continue;
		VLW vlw1;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str())-1;
		vlw1.chrEnd=atoi(vs[2].c_str());
		vlw1.null=true;
		vlw1.Rabl=atof(vs[3].c_str())-0.5;
		vlw1.blackListed=false;
		vlw1.Excluded=-1;
//		vlw1.Rabl_index=nRabl*vlw1.Rabl;
		vlws.push_back(vlw1);
	}
}

int findInVLWs(const VLW& vlw1, const vector<VLW>& vlws, int start, int end)
{
	if (start>=end) return -1;
	int j=(start+end)/2;
	if (vlw1.chrName==vlws[j].chrName && vlw1.chrStart>=vlws[j].chrStart && vlw1.chrEnd<=vlws[j].chrEnd) return j;
	else if (leftNonOverlap(vlw1,vlws[j])) return findInVLWs(vlw1,vlws,start,j);
	else if (leftNonOverlap(vlws[j],vlw1)) return findInVLWs(vlw1,vlws,j+1,end);
	else return -1;
}

int vlwdis(const vector<VLW>& vlws, int i, int j)
{
	if (vlws[i].chrName!=vlws[j].chrName) return -1;
	if (i>j) swap(i,j);
	return (vlws[j].chrStart+vlws[j].chrEnd-vlws[i].chrStart-vlws[i].chrEnd)/2;
}

void AssignVLWs(vector<VLW>& vlws, const vector<VLW>& vlws1) //assign color and VLW number for each bin in vlws
{
	for (int i=0;i<vlws.size();++i)
	{
		int j=findInVLWs(vlws[i],vlws1,0,vlws1.size());
		vlws[i].vlwsnumber=j;
		if (j>=0) vlws[i].color=vlws1[j].color;
		else vlws[i].color=-1;
	}
}
int findInvlws(int j,const vector<int>& idxvlws,int start,int end)
{
	if (start>=end) return -1;
	int k=(start+end)/2;
	if (j==idxvlws[k]) return k;
	else if (j>idxvlws[k]) return findInvlws(j,idxvlws,k+1,end);
	else return findInvlws(j,idxvlws,start,k);
}

class Analysis
{
public:
	double d0;
	int dmin;
	string chrana;
//	double steplength;
	int anaTp;
	int nsession;
	vector<VLW> regions;
	map<string,pair<int,int> > chrIdxRange;
	vector<map<int,int> > interactions;
	vector<vector<Block> > excludedBlocks;
	vector<double> bias;
	vector<vector<double> > cscore;
	double ralfa, rbeta, rgamma; //rabl parameters
//	vector<vector<double> > rabl; //rabl factors
	double maxdevc;
	vector<double> hh;
//	vector<double> h;
	vector<double> Nd;
	vector<double> Nr;
//	vector<vector<double> > Nrr; //Number of interactions between two Rabl index categories
	double Nt;
	double delta;
//	map<string,double> sumchr;
	
//	double sumb;
	Matrix sumc; //compartment*3, c, c*Rabl, c*Rabl^2 
	map<string,Matrix > sumcchr; //chromosome, compartment*3
	vector<vector<Matrix> > sumcBlocks;//Blockgroup, Block, compartment*3
//	map<string,vector<vector<double> > > sumcchrrabl;//string: chrName, sumcchrrabl[chrName][k][l]: k=0,...,Nsub-1; l=0,...,nRabl-1. sum of all Cscore[k] for windows with Rabl_index l on chr chrName.
	double likelihood;
	int readRegions(const string& fileName);
	int readBlackList(const string& fileName);
	int readExcludedBlocks(const string& fileName);
	int readInteractions(const string& fileName);
	int readMatrix(const string& fileName);
	int initBias();
	int initCscore();
	int updatehh();
	int updateBias();
	int updateCscore();
	int updateRabl();
	int updateLikelihood();
	int normalize();
};

int Analysis::readRegions(const string& fileName)
{
	readVLWindow(fileName,regions);
	cout<<regions.size()<<endl;
	cout<<"read VLWs File"<<endl;
	sort(regions.begin(),regions.end(),leftNonOverlap);
	cout<<"VLWs sorted"<<endl;
	d0=regions[0].chrEnd-regions[0].chrStart;
	chrIdxRange.clear();
	for (int i=0;i<regions.size();++i)
	{
		string cn=regions[i].chrName;
		map<string,pair<int,int> >::iterator itr1=chrIdxRange.find(cn);
		if (itr1==chrIdxRange.end()) 
		{
			chrIdxRange.insert(make_pair(cn,make_pair(i,i+1)));
		}
		else itr1->second.second=i+1;
	}
	cout<<"chr ranges assigned"<<endl;
	excludedBlocks.clear();
	return 0;
}

int Analysis::readBlackList(const string& fileName)
{
	ifstream fin(fileName.c_str());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<3) continue;
		Block blk1;
		blk1.chrName=vs[0];
		blk1.chrStart=atoi(vs[1].c_str());
		blk1.chrEnd=atoi(vs[2].c_str());
		int start1=chrIdxRange[blk1.chrName].first;
		int end1=chrIdxRange[blk1.chrName].second;
		blk1.blockStart=start1;
		blk1.blockEnd=end1;
		for (int i=start1;i<end1;++i)
		{
			if (regions[i].chrEnd>blk1.chrStart)
			{
				blk1.blockStart=i;
				break;	
			}	
		}
		for (int i=end1-1;i>=start1;--i)
		{
			if (regions[i].chrStart<blk1.chrEnd)
			{
				blk1.blockEnd=i+1;
				break;	
			}	
		}
		for (int j=blk1.blockStart;j<blk1.blockEnd;++j)
		{
			regions[j].blackListed=true;
		}
	}
}

int Analysis::readExcludedBlocks(const string& fileName)
{
	ifstream fin(fileName.c_str());
	excludedBlocks.clear();
	int ct=0;
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs1;
		parseStr(s,'\t',vs1);
		if (vs1.size()<2) continue;
		vector<Block> vb(0);
		for (int i=0;i<vs1.size();++i)
		{
			vector<string> vs;
			string s1=vs1[i];
			parseStr(s1,'_',vs);
			Block blk1;
			blk1.chrName=vs[0];
			blk1.chrStart=atoi(vs[1].c_str());
			blk1.chrEnd=atoi(vs[2].c_str());
			int start1=chrIdxRange[blk1.chrName].first;
			int end1=chrIdxRange[blk1.chrName].second;
			blk1.blockStart=start1;
			blk1.blockEnd=end1;
			for (int i=start1;i<end1;++i)
			{
				if (regions[i].chrEnd>blk1.chrStart)
				{
					blk1.blockStart=i;
					break;	
				}	
			}
			for (int i=end1-1;i>=start1;--i)
			{
				if (regions[i].chrStart<blk1.chrEnd)
				{
					blk1.blockEnd=i+1;
					break;	
				}	
			}
			vb.push_back(blk1);
			for (int j=blk1.blockStart;j<blk1.blockEnd;++j)
			{
				regions[j].Excluded=ct;
			}
		}
		excludedBlocks.push_back(vb);
		++ct;
	}
}

int Analysis::readInteractions(const string& fileName)
{
	ifstream fin(fileName.c_str());
	if (!fin.is_open())
	{
		cout<<"Bad interaction File Name"<<endl;
		return -1;	
	}
	int count=0;
	interactions.clear();
	Nt=0;
	Nd.clear();
	Nd.resize(132,0);
	hh.clear();
	hh.resize(132,0);
	Nr.clear();
	Nr.resize(regions.size(),0);
//	Nrr.clear();
//	Nrr.resize(nRabl,vector<double>(nRabl,0));
	interactions.clear();
	interactions.resize(regions.size());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<6) continue;
		++count;
		if (count%1000000==0) cout<<count<<endl;
		VLW vlw1;
		VLW vlw2;
		vlw1.chrName=vs[1];
		vlw1.chrStart=atoi(vs[2].c_str());
		vlw1.chrEnd=vlw1.chrStart+1;
		vlw2.chrName=vs[4];
		vlw2.chrStart=atoi(vs[5].c_str());
		vlw2.chrEnd=vlw2.chrStart+1;	
		int i1=findInVLWs(vlw1,regions,0,regions.size());
		int i2=findInVLWs(vlw2,regions,0,regions.size());
		if (i1<0 ||i2<0) continue;
		if (regions[i1].blackListed || regions[i2].blackListed) continue;
		if (regions[i1].chrName==regions[i2].chrName) continue; // Only use inter-chr reads
			if (regions[i1].Excluded>=0 && regions[i1].Excluded==regions[i2].Excluded) continue; //both i1 and i2 in same excluded blocks, don't count it.
		map<int,int>::iterator itr1=interactions[i1].find(i2);
		if (itr1==interactions[i1].end()) interactions[i1].insert(make_pair(i2,1));
		else itr1->second++;
		map<int,int>::iterator itr2=interactions[i2].find(i1);
		if (itr2==interactions[i2].end()) interactions[i2].insert(make_pair(i1,1));
		else itr2->second++;
		regions[i1].null=false;
		regions[i2].null=false;
		Nr[i1]++;
		Nr[i2]++;
		++Nt;
//		++Nrr[regions[i1].Rabl_index][regions[i2].Rabl_index];
//		++Nrr[regions[i2].Rabl_index][regions[i1].Rabl_index];
	}
	/*ofstream flog("interactions.txt");
	for (int i=0;i<interactions.size();++i)
	{
		for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
		{
			if (itr->first>i) break;
			flog<<i<<"\t"<<itr->first<<"\t"<<itr->second<<endl;
		}	
	}*/
	return 0;
}

/*//This function is to make sumcc==total reads
int Analysis::normalize()
{
	double sum1=0,sum2=0;
	for (int u=0;u<nRabl;++u)
	{
		for (int v=0;v<nRabl;++v)
		{
			for (int k=0;k<Nsub;++k)
			{
				sum1+=sumc[k][u]*sumc[k][v];
				for (map<string,vector<vector<double> > >::iterator itr=sumcchr.begin();itr!=sumcchr.end();++itr)
				{
					sum1-=itr->second[k][u]*itr->second[k][v];	
				}	
			}	
			sum2+=Nrr[u][v];
		}	
	}
	double factor=sqrt(sum2/sum1);
	for (int i=0;i<cscore.size();++i)
	{
		for (int k=0;k<Nsub;++k)
		{
			cscore[i][k]*=factor;	
		}
	}
	for (int u=0;u<nRabl;++u)
	{
		for (int v=0;v<nRabl;++v)
		{
			rabl[u][v]*=sum1/sum2;	
		}	
	}
	return 0;
}*/

/*
int Analysis::readMatrix(const string& fileName)
{
	ifstream fin(fileName.c_str());
	int count=0;
	interactions.clear();
	Nd.clear();
	Nd.resize(132,0);
	hh.clear();
	hh.resize(132,0);
	Nr.clear();
	Nr.resize(regions.size(),0);
	interactions.clear();
	interactions.resize(regions.size());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<4) continue;
		++count;
		if (count%1000000==0) cout<<count<<endl;
		VLW vlw1;
		VLW vlw2;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str());
		vlw1.chrEnd=vlw1.chrStart+1;
		vlw2.chrName=vs[0];
		vlw2.chrStart=atoi(vs[2].c_str());
		vlw2.chrEnd=vlw2.chrStart+1;
		double x=atof(vs[3].c_str());
		int i1=findInVLWs(vlw1,regions,0,regions.size());
		int i2=findInVLWs(vlw2,regions,0,regions.size());
		if (i1<0 ||i2<0) continue;
		if (regions[i1].chrName!=regions[i2].chrName) continue;
		if (chrana!="" && (regions[i1].chrName!=chrana ||regions[i2].chrName!=chrana)) continue;
		double d=d0*abs(i2-i1);
		int dd=(log10(d+0.01)-3)/steplength;
		if (dd<dmin) continue;
		if (dd>=hh.size()) dd=hh.size()-1;
		map<int,int>::iterator itr1=interactions[i1].find(i2);
		if (itr1==interactions[i1].end()) interactions[i1].insert(make_pair(i2,1));
		else itr1->second+=x;
		map<int,int>::iterator itr2=interactions[i2].find(i1);
		if (itr2==interactions[i2].end()) interactions[i2].insert(make_pair(i1,1));
		else itr2->second+=x;
		Nd[dd]++;
		Nr[i1]++;
		Nr[i2]++;
	}
	return 0;
}*/
/*
int Analysis::initBias()
{
	bias.clear();
	bias.resize(regions.size(),0);
	for (int i=0;i<bias.size();++i)
	{
		bias[i]=Nr[i]/sqrt(2.0*Nt);	
	}
	sumchr.clear();
	sumb=0;
	for (int i=0;i<regions.size();++i)
	{
		sumb+=bias[i];
		map<string,double>::iterator itr=sumchr.find(regions[i].chrName);
		if (itr==sumchr.end())
		{
			sumchr.insert(make_pair(regions[i].chrName,bias[i]));	
		}
		else itr->second+=bias[i];
	}
	return 0;	
}*/

int Analysis::initCscore()
{
	cscore.clear();
	cscore.resize(regions.size(),vector<double>(Nsub,0));
	ralfa=0;
	rbeta=0;
//	rgamma=0;
	for (int i=0;i<cscore.size();++i)
	{
		if (interactions[i].size()>=thr_interaction) raninit(cscore[i]);
	}
	maxdevc=0;
	sumc.clear();
	sumc.resize(Nsub,vector<double>(3,0));
	sumcchr.clear();
	for (int i=0;i<regions.size();++i)
	{
		double r=regions[i].Rabl;
		for (int k=0;k<Nsub;++k)
		{
			sumc[k][0]+=cscore[i][k];
			sumc[k][1]+=cscore[i][k]*r;
			sumc[k][2]+=cscore[i][k]*r*r;
		}
		map<string,Matrix>::iterator itr=sumcchr.find(regions[i].chrName);
		if (itr==sumcchr.end())
		{
			Matrix vv(Nsub,vector<double>(3,0));
			for (int k=0;k<Nsub;++k)
			{
				vv[k][0]=cscore[i][k];
				vv[k][1]=cscore[i][k]*r;
				vv[k][2]=cscore[i][k]*r*r;	
			}
			sumcchr.insert(make_pair(regions[i].chrName,vv));	
		}
		else 
		{
			for (int k=0;k<Nsub;++k)
			{
				itr->second[k][0]+=cscore[i][k];
				itr->second[k][1]+=cscore[i][k]*r;
				itr->second[k][2]+=cscore[i][k]*r*r;
			}	
		}
	}
	sumcBlocks.clear();
	for (int l=0;l<excludedBlocks.size();++l)
	{
		vector<Matrix> sum1(excludedBlocks[l].size(),vector<vector<double> >(Nsub,vector<double>(3,0)));
		for (int i=0;i<sum1.size();++i)
		{
			for (int j=excludedBlocks[l][i].blockStart;j<excludedBlocks[l][i].blockEnd;++j)
			{
				for (int k=0;k<Nsub;++k)
				{
					sum1[i][k][0]+=cscore[j][k];
					sum1[i][k][1]+=cscore[j][k]*regions[j].Rabl;
					sum1[i][k][2]+=cscore[j][k]*regions[j].Rabl*regions[j].Rabl;
				}
			}
		}
		sumcBlocks.push_back(sum1);
	}
	return 0;		
}
/*
int Analysis::updatehh()
{
	for ()
//	hh[dmin]=1;
	for (int i=dmin;i<hh.size();++i)
	{
		if (Nd[i]>0) hh[i]=Nd[i]/sumcc[i];
	}
	sumbh.clear();
	sumbh.resize(regions.size(),0);
	sumch.clear();
	sumch.resize(regions.size(),0);
	sumch1.clear();
	sumch1.resize(regions.size(),0);
	#pragma omp parallel for
	for (int i=0;i<regions.size();++i)
	{
		if (interactions[i].size()<=1) continue;
		string cn=regions[i].chrName;
		int r1=chrIdxRange[cn].first;
		int r2=chrIdxRange[cn].second;
		double F=0,G=0,G1=0;
		for (int j=r1;j<r2;++j)
		{
			double d=abs(j-i)*d0;
			int dd=(log10((double)d+0.01)-3)/steplength;
			if (dd<dmin) continue;
			if (dd>=hh.size()) dd=hh.size()-1;
			double u=bias[j]*hh[dd];
			F+=u;
			G+=u*cscore[j];
			G1+=u*cscore1[j];
		}
		sumbh[i]=F;
		sumch[i]=G;
		sumch1[i]=G1;
	}
	double L=0;
	for (int i=0;i<regions.size();++i)
	{
		if (interactions[i].size()<=1) continue;
		L+=Nr[i]*log(bias[i]);
		for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
		{
			int j=itr->first;
			if (j>i) break;
			double nij=itr->second;
			L+=nij*log(1+cscore[i]*cscore[j]+cscore1[i]*cscore1[j]);
		}
	}
	for (int i=dmin;i<hh.size();++i)
	{
		if (hh[i]>0) L+=Nd[i]*log(hh[i]);
	}
//	L-=sumcc[dmin];
	cout<<"L="<<L<<"\tmaxdevc="<<maxdevc<<endl;
	if (L<likelihood+1) return 1;
	likelihood=L;
	return 0;	
}*/

/*
int Analysis::updateLikelihood()
{
	double L=0;
//	int sum1=0;
	for (int i=0;i<regions.size();++i)
	{
//		if (interactions[i].size()<thr_interaction) continue;
		for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
		{
			int j=itr->first;
			if (j>i) break;
			double nij=itr->second;
			double cc=0;
			for (int k=0;k<Nsub;++k)
			{
				cc+=cscore[i][k]*cscore[j][k];
			}
			L+=nij*log(cc);
//			sum1+=nij;
		}
	}
//	cout<<"sum1="<<sum1<<endl;

	vector<vector<double> > sumcc(nRabl,vector<double>(nRabl,0));
	
//	int sum2=0;
	for (int u=0;u<nRabl;++u)
	{
		for (int v=0;v<nRabl;++v)
		{
			for (int k=0;k<Nsub;++k)
			{
				sumcc[u][v]+=sumc[k][u]*sumc[k][v];
				for (map<string,vector<vector<double> > >::iterator itr=sumcchr.begin();itr!=sumcchr.end();++itr)
				{
					sumcc[u][v]-=itr->second[k][u]*itr->second[k][v];	
				}
			}
			rabl[u][v]=Nrr[u][v]/sumcc[u][v];
			if (rabl[u][v]>0) L+=(Nrr[u][v]*log(rabl[u][v])-Nrr[u][v])/2;
//			sum2+=Nrr[u][v];	
			cout<<rabl[u][v]<<"\t";
		}
		cout<<endl;
	}	
//	cout<<"sum2="<<sum2<<endl;
	
	cout<<"L="<<L<<endl;
	if (L<likelihood+0.1) return 1;
	likelihood=L;
	return 0;
}*/
/*
int Analysis::updateBias()
{
	for (map<string,pair<int,int> >::iterator itr=chrIdxRange.begin();itr!=chrIdxRange.end();++itr)
	{
		string cn=itr->first;
		for (int i=itr->second.first;i<itr->second.second;++i)
		{
		//	if (interactions[i].size()<thr_interaction) continue;
			double bf=0;//sumb-sumchr[cn];
			for (int k=0;k<Nsub;++k)
			{
				bf+=cscore[i][k]*(sumc[k]-sumcchr[cn][k]);
			}
			bias[i]=Nr[i]/bf;
		}
		sumchr[cn]=0;
		sumcchr[cn].clear();
		sumcchr[cn].resize(Nsub,0);
		for (int i=itr->second.first;i<itr->second.second;++i)
		{
			sumchr[cn]+=bias[i];
			for (int k=0;k<Nsub;++k)
			{
				sumcchr[cn][k]+=bias[i]*cscore[i][k];	
			}
		}
		sumb=0;
		for (map<string,double>::iterator itr1=sumchr.begin();itr1!=sumchr.end();++itr1)
		{
			sumb+=itr1->second;
		}
		for (int k=0;k<Nsub;++k)
		{
			sumc[k]=0;	
		}
		for (map<string,vector<double> >::iterator itr2=sumcchr.begin();itr2!=sumcchr.end();++itr2)
		{
			for (int k=0;k<Nsub;++k)
			{
				sumc[k]+=itr2->second[k];
			}
		}
	}
	//cout<<"UpdateBiasFinished"<<endl;
}*/

int Analysis::updateLikelihood()
{
	double A=0,B=0,C=0,D=0;
	for (map<string,Matrix>::iterator itr=sumcchr.begin();itr!=sumcchr.end();++itr)
	{
		for (int k=0;k<Nsub;++k)
		{
			A+=itr->second[k][0]*(sumc[k][0]-itr->second[k][0])/2;
//			B+=itr->second[k][1]*(sumc[k][0]-itr->second[k][0]);
//			C+=itr->second[k][2]*(sumc[k][0]-itr->second[k][0]);
			B+=itr->second[k][1]*(sumc[k][1]-itr->second[k][1])/2;
			C+=itr->second[k][2]*(sumc[k][1]-itr->second[k][1]);
		}
	}
	for (int l=0;l<sumcBlocks.size();++l)
	{
		for (int k=0;k<Nsub;++k)
		{
			for (int m1=0;m1<sumcBlocks[l].size();++m1)
			{
				for (int m2=m1+1;m2<sumcBlocks[l].size();++m2)
				{
					if (excludedBlocks[l][m1].chrName==excludedBlocks[l][m2].chrName) continue;
					A-=sumcBlocks[l][m1][k][0]*sumcBlocks[l][m2][k][0];
					B-=sumcBlocks[l][m1][k][1]*sumcBlocks[l][m2][k][1];
					C-=sumcBlocks[l][m1][k][1]*sumcBlocks[l][m2][k][2]+sumcBlocks[l][m1][k][2]*sumcBlocks[l][m2][k][1];
				}
			}	
		}
		
	}
//	cout<<"A="<<A<<"\tB="<<B<<"\tC="<<C<<"\tD="<<D<<endl;
	vector<double> x(2,0);
	x[0]=ralfa;
	x[1]=rbeta;
	double eps=1.0e-8;
	//Newton search
	for (int nlp1=0;nlp1<100;++nlp1)
	{
		vector<double> FL(2,0);
		Matrix GL(2,vector<double>(2,0));
		for (int i=0;i<interactions.size();++i)
		{
			for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
			{
				int j=itr->first;
				if (j>i) break;
				register int nij=itr->second;
				register double u=regions[i].Rabl*regions[j].Rabl;
				register double v=regions[i].Rabl+regions[j].Rabl;
				register double Rij=1+u*(x[0]+x[1]*v);
				register double u_Rij=u/Rij;
				register double uv_Rij=u/Rij*v;
				register double nu_Rij=nij*u_Rij;
				register double nuv_Rij=nij*uv_Rij;
				FL[0]+=nu_Rij;
				FL[1]+=nuv_Rij;
				GL[0][0]+=nu_Rij*u_Rij;
				GL[0][1]+=nu_Rij*uv_Rij;
				GL[1][1]+=nuv_Rij*uv_Rij;
			}
		}
		GL[1][0]=GL[0][1];
		FL[0]-=B;
		FL[1]-=C;
		vector<double> dx=invp(GL,FL);
		double rd=sqrt(innerprod(FL,dx)/2); //gradient
//		cout<<x[0]<<"\t"<<x[1]<<"\t"<<dx[0]<<"\t"<<dx[1]<<"\trd="<<rd<<endl;
		if (rd<eps) break; //stops
		else
		{
			x[0]+=dx[0];
			x[1]+=dx[1];
		}
	}
	ralfa=x[0];
	rbeta=x[1];
	
	double L=0;
	for (int i=0;i<regions.size();++i)
	{		for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
		{
			int j=itr->first;
			if (j>i) break;
			double nij=itr->second;
			double cc=0;
			for (int k=0;k<Nsub;++k)
			{
				cc+=cscore[i][k]*cscore[j][k];
			}
			double u=regions[i].Rabl*regions[j].Rabl;
			double v=regions[i].Rabl+regions[j].Rabl;
			double Rij=1+u*(ralfa+rbeta*v);
			L+=nij*log(cc*Rij);
		}
	}
	L-=A+ralfa*B+rbeta*C;
	cout<<"L="<<L<<"\tralfa="<<ralfa<<"\trbeta="<<rbeta<<endl;
	if (L<likelihood+1) return 1;
	likelihood=L;
	return 0;
}
/*
int Analysis::updateBias()
{
	for (map<string,pair<int,int> >::iterator itr=chrIdxRange.begin();itr!=chrIdxRange.end();++itr)
	{
		string cn=itr->first;
		for (int i=itr->second.first;i<itr->second.second;++i)
		{
		//	if (interactions[i].size()<thr_interaction) continue;
			double bf=0;//sumb-sumchr[cn];
			for (int k=0;k<Nsub;++k)
			{
				bf+=cscore[i][k]*(sumc[k]-sumcchr[cn][k]);
			}
			bias[i]=Nr[i]/bf;
		}
		sumchr[cn]=0;
		sumcchr[cn].clear();
		sumcchr[cn].resize(Nsub,0);
		for (int i=itr->second.first;i<itr->second.second;++i)
		{
			sumchr[cn]+=bias[i];
			for (int k=0;k<Nsub;++k)
			{
				sumcchr[cn][k]+=bias[i]*cscore[i][k];	
			}
		}
		sumb=0;
		for (map<string,double>::iterator itr1=sumchr.begin();itr1!=sumchr.end();++itr1)
		{
			sumb+=itr1->second;
		}
		for (int k=0;k<Nsub;++k)
		{
			sumc[k]=0;	
		}
		for (map<string,vector<double> >::iterator itr2=sumcchr.begin();itr2!=sumcchr.end();++itr2)
		{
			for (int k=0;k<Nsub;++k)
			{
				sumc[k]+=itr2->second[k];
			}
		}
	}
	//cout<<"UpdateBiasFinished"<<endl;
}*/

int Analysis::updateCscore()
{	
	maxdevc=0;
	for (map<string,pair<int,int> >::iterator itr=chrIdxRange.begin();itr!=chrIdxRange.end();++itr)
	{
		string cn=itr->first;
//		cout<<cn<<endl;
		#pragma omp parallel for
		for (int i=itr->second.first;i<itr->second.second;++i)
		{
			if (interactions[i].size()<thr_interaction) continue;
			vector<double> F(Nsub,0);
			for (int k=0;k<Nsub;++k)
			{
				double r=regions[i].Rabl;
				F[k]=sumc[k][0]-sumcchr[cn][k][0]+(ralfa+rbeta*r)*r*(sumc[k][1]-sumcchr[cn][k][1])+rbeta*r*(sumc[k][2]-sumcchr[cn][k][2]);
				int l=regions[i].Excluded;
				if (l>=0)
				{
					for (int m=0;m<sumcBlocks[l].size();++m)
					{
						if (excludedBlocks[l][m].chrName==cn) continue;
						F[k]-=sumcBlocks[l][m][k][0]+(ralfa+rbeta*r)*r*sumcBlocks[l][m][k][1]+rbeta*r*sumcBlocks[l][m][k][2];
					}	
				}
			}
			vector<double> x=cscore[i];
			vector<double> lb(Nsub,1.0);
			double t=1.0,yeta=30,mu=10,eps=1.0e-8;
			
			//get initial Likelihood associated...
			double L0=0;
			for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
			{
				int j=itr->first;
				int nij=itr->second;
				double sxc=0;
				for (int k=0;k<Nsub;++k)
				{
					sxc+=x[k]*cscore[j][k];
				}
				L0+=nij*log(sxc);
			}
			for (int k=0;k<Nsub;++k)
			{
				L0-=F[k]*x[k];
			}
//			cout<<"L0="<<L0<<endl;
		
//			cout<<"i="<<i<<"\t"<<interactions[i].size()<<endl;
			for (int nlp1=0;nlp1<100;++nlp1)
			{
//				cout<<"nlp1="<<nlp1<<endl;
//				cout<<x[0]<<"\t"<<x[1]<<"\t"<<lb[0]<<"\t"<<lb[1]<<"\t"<<nu<<endl;
				Matrix GL(Nsub,vector<double>(Nsub,0));
				vector<double> RX(Nsub,0);
				vector<double> RL(Nsub,0);
				//get Newton search direction
				int j,nij;
				vector<double> u(Nsub,0);
				vector<double> y(Nsub,0);
				vector<vector<double> > z(Nsub,vector<double>(Nsub,0));
				for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
				{
					j=itr->first;
					nij=itr->second;
					double sxc=0;
					for (int k=0;k<Nsub;++k)
					{
						sxc+=x[k]*cscore[j][k];
					}
					for (int k=0;k<Nsub;++k)
					{
						u[k]=cscore[j][k]/sxc;
						y[k]=nij*u[k];
						RX[k]+=y[k];	
					}
					for (int k=0;k<Nsub;++k)
					{
						for (int l=0;l<Nsub;++l)
						{
							z[k][l]=y[k]*u[l];
							GL[k][l]+=z[k][l];
						}
					}
				}
//				cout<<"GL="<<"\t"<<GL[0][0]<<"\t"<<GL[0][1]<<"\t"<<GL[1][0]<<"\t"<<GL[1][1]<<endl;
				for (int k=0;k<Nsub;++k)
				{
					RX[k]+=lb[k]-F[k];
					RL[k]=1/t-lb[k]*x[k];
					GL[k][k]+=lb[k]/x[k];
				}
				vector<double> FL(RX);
				for (int k=0;k<Nsub;++k)
				{
					FL[k]+=RL[k]/x[k];
				}
				vector<double> dx=invp(GL,FL);
				vector<double> dl(Nsub,0);
				for (int k=0;k<Nsub;++k)
				{
					dl[k]=(RL[k]-lb[k]*dx[k])/x[k];
				}
			/*	for (int k=0;k<dl.size();++k)
				{
					cout<<dl[k]<<"\t";
				}
				cout<<endl;*/
//				cout<<dx[0]<<"\t"<<dx[1]<<"\t"<<dl[0]<<"\t"<<dl[1]<<"\t"<<dn<<endl;
				//back track search to get the step length
				double s=0.9999,alpha=0.05,beta=0.5;
				for (;;) //feasibility
				{
					vector<double> xx(Nsub,0);
					vector<double> ll(Nsub,0);
					int flag=0;
					for (int k=0;k<Nsub;++k)
					{
						xx[k]=x[k]+s*dx[k];
						ll[k]=lb[k]+s*dl[k];
						if (xx[k]<=0 || ll[k]<=0)
						{
							flag=1;
							break;
						}
					}
					if (flag>0) s*=beta;
					else break;
				}
//				cout<<"s="<<s<<endl;
				for (int nsp=0;nsp<100;++nsp)
				{
				//	cout<<"nsp="<<nsp<<endl;
					vector<double> RL1(Nsub,0);
					vector<double> RX1(Nsub,0);
					vector<double> xx(Nsub,0);
					vector<double> ll(Nsub,0);
					int flag=0;
					for (int k=0;k<Nsub;++k)
					{
						xx[k]=x[k]+s*dx[k];
						ll[k]=lb[k]+s*dl[k];
					}
					for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
					{
						j=itr->first;
						nij=itr->second;
						double sxc=0;
						for (int m=0;m<Nsub;++m)
						{
							sxc+=xx[m]*cscore[j][m];
						}
						for (int k=0;k<Nsub;++k)
						{
							u[k]=cscore[j][k]/sxc;
							y[k]=nij*u[k];
							RX1[k]+=y[k];	
						}
					}
					for (int k=0;k<Nsub;++k)
					{
						RX1[k]+=ll[k]-F[k];
						RL1[k]=-ll[k]*xx[k]+1/t;
					}
					double aa=1-alpha*s;
//					cout<<"R="<<RX[0]<<"\t"<<RX[1]<<"\t"<<RL[0]<<"\t"<<RL[1]<<"\t"<<RN<<endl;
//					cout<<"R1="<<RX1[0]<<"\t"<<RX1[1]<<"\t"<<RL1[0]<<"\t"<<RL1[1]<<"\t"<<RN1<<endl;
					if ((innerprod(RX1)+innerprod(RL1))<=aa*aa*(innerprod(RX)+innerprod(RL)))
					{
						x=xx;
						lb=ll;
						RX=RX1;
						RL=RL1;
						break;	
					}
					else s*=beta;
//					cout<<xx[0]<<"\t"<<xx[1]<<"\t"<<ll[0]<<"\t"<<ll[1]<<"\t"<<nn<<endl;
				}
				//back track search end, see the result
				double rd=sqrt(innerprod(RX,RX)+innerprod(RL,RL)); //gradient
//				cout<<"rd="<<rd<<endl;
				yeta=innerprod(x,lb);
				if (rd<eps && yeta<eps) break; //stops
				else //update t, redo the search
				{
					t=4.0*mu/yeta;
			//		cout<<t<<endl;
				}
			}
			
			double L1=0;
			for (map<int,int>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
			{
				int j=itr->first;
				int nij=itr->second;
				double sxc=0;
				for (int k=0;k<Nsub;++k)
				{
					sxc+=x[k]*cscore[j][k];
				}
				L1+=nij*log(sxc);
			}
			for (int k=0;k<Nsub;++k)
			{
				L1-=F[k]*x[k];
			}
			if (L1<L0) cout<<cn<<"\t"<<i<<"\tL0="<<L0<<"\tL1="<<L1<<endl;
			cscore[i]=x;
		}
		for (int k=0;k<Nsub;++k)
		{
			sumcchr[cn][k][0]=0;
			sumcchr[cn][k][1]=0;
			sumcchr[cn][k][2]=0;

			for (int i=itr->second.first;i<itr->second.second;++i)
			{
				double r=regions[i].Rabl;
				sumcchr[cn][k][0]+=cscore[i][k];
				sumcchr[cn][k][1]+=cscore[i][k]*r;
				sumcchr[cn][k][2]+=cscore[i][k]*r*r;
			}
		}
		for (int k=0;k<Nsub;++k)
		{
			sumc[k][0]=0;
			sumc[k][1]=0;
			sumc[k][2]=0;
			for (map<string,Matrix>::iterator itr2=sumcchr.begin();itr2!=sumcchr.end();++itr2)
			{
				sumc[k][0]+=itr2->second[k][0];
				sumc[k][1]+=itr2->second[k][1];
				sumc[k][2]+=itr2->second[k][2];
			}	
		}
		for (int l=0;l<excludedBlocks.size();++l)
		{
			
			for (int m=0;m<excludedBlocks[l].size();++m)
			{
				if (excludedBlocks[l][m].chrName!=cn) continue;
				Matrix sum1(Nsub,vector<double>(3,0));
				for (int j=excludedBlocks[l][m].blockStart;j<excludedBlocks[l][m].blockEnd;++j)
				{
					for (int k=0;k<Nsub;++k)
					{
						sum1[k][0]+=cscore[j][k];
						sum1[k][1]+=cscore[j][k]*regions[j].Rabl;
						sum1[k][2]+=cscore[j][k]*regions[j].Rabl*regions[j].Rabl;
					}
				}
				sumcBlocks[l][m]=sum1;
			}
			
		}
	}
//	cout<<"UpdateCscoreFinished"<<endl;
}

int main(int argc, char* argv[])
{
	if (argc<7)
	{
		cout<<"Usage: CscoreTool-M <windows_Rabl.bed> <input.summary> <OutputPrefix> <N_subcompartments> <N_Rablwindow> <session> [Blacklist.bed] [ExcludedInteractions.txt]"<<endl;
		return 0;	
	}
	cout<<setprecision(17);
	SetSeed(CreateSeed());
	string chrs="";
	Analysis analysis;
	analysis.chrana=chrs;
	nRabl=atoi(argv[5]);
	analysis.nsession=atoi(argv[6]);

	Nsub=atoi(argv[4]);
	cout<<"Nsub="<<Nsub<<endl;
//	int minDis=atoi(argv[5]);
//	analysis.dmin=(log10((double)minDis+0.01)-3)/steplength;
//	cout<<"dmin="<<analysis.dmin<<endl;

	string outputprefix=argv[3];
	string windowFileName=argv[1];
	int f1=analysis.readRegions(windowFileName);
	if (f1<0) return 0;
	if (argc>7)
	{
		string blackListFileName=argv[7];
		analysis.readBlackList(blackListFileName);	
	}	
		
	if (argc>8)
	{
		string excludedBlockFileName=argv[8];
		analysis.readExcludedBlocks(excludedBlockFileName);
	}
	string inputFileName=argv[2];
	int f2=analysis.readInteractions(inputFileName);
	if (f2<0) return 0;
	int count=0;
//	analysis.initBias();
	analysis.initCscore();
	analysis.likelihood=-1.0e308;
	omp_set_num_threads(analysis.nsession);
	for (int i=0;i<500;++i)
	{
		analysis.updateCscore();	
		int flg=analysis.updateLikelihood();
		if (flg>0) break;		
//		analysis.updateBias();
		
	}
	
	
	analysis.bias.resize(analysis.cscore.size(),0);
	for (int i=0;i<analysis.bias.size();++i)
	{
		analysis.bias[i]=0;
		for (int k=0;k<Nsub;++k)
		{
			analysis.bias[i]+=analysis.cscore[i][k];	
		}
	}
	
	string ofname0=outputprefix+"bias.txt";
	ofstream fbias(ofname0.c_str());
	for (int i=0;i<analysis.bias.size();++i)
	{
		fbias<<i;//idxvlws[i];
		fbias<<"\t"<<analysis.bias[i];
		fbias<<endl;
	}
	
	string ofname1=outputprefix+"_cscore.txt";
	ofstream fegc(ofname1.c_str());
	for (int i=0;i<analysis.cscore.size();++i)
	{
		fegc<<i;//idxvlws[i];
		for (int k=0;k<Nsub;++k)
		{
			fegc<<"\t"<<analysis.cscore[i][k];
		}
		fegc<<endl;
	}
	
	string ofname2=outputprefix+"_cscore.bedgraph";
	ofstream fbdg(ofname2.c_str());
	for (int k=0;k<Nsub;++k)
	{
		string trackname=outputprefix+"_"+char('0'+k);
		fbdg<<"track type=\"bedGraph\" name=\""<<trackname<<"\" visibility=full viewLimits=0:1 windowingFunction=mean autoScale=off"<<endl;
		for (int i=0;i<analysis.cscore.size();++i)
		{
				fbdg<<analysis.regions[i].chrName<<"\t"<<analysis.regions[i].chrStart<<"\t"<<analysis.regions[i].chrEnd;
				if (analysis.bias[i]>0) fbdg<<"\t"<<analysis.cscore[i][k]/analysis.bias[i];
				else fbdg<<"\t"<<0;
				fbdg<<endl;
		}
	}
	
	string ofname3=outputprefix+"_rabl.txt";
	ofstream frabl(ofname3.c_str());
	frabl<<"# "<<analysis.ralfa<<"\t"<<analysis.rbeta<<"\t"<<analysis.rgamma<<endl;
	for (double u=0.5;u<nRabl;++u)
	{
		for (double v=0.5;v<nRabl;++v)
		{
			double x=u/nRabl-0.5;
			double y=v/nRabl-0.5;
			double r=1+(analysis.ralfa+analysis.rbeta*(x+y))*x*y;
			frabl<<u/nRabl<<"\t"<<v/nRabl<<"\t"<<r<<endl;
		}
	}
}
