
using namespace std;
#include <math.h>
//#include <map>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
//#include <algorithm>
#include <sstream>

//#include <gsl/gsl_permutation.h>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

struct assignment{
	int locus;
	int chr;
	string name;
	vector <double> prob;
};


void InputData(vector <assignment> &assignmentV, std::string infile1);


int main(int argc, char *argv[]){

	string infile1(argv[1]); // assignment prob file with a two lines header*/
	
	vector <assignment> assignmentV;
	InputData(assignmentV,infile1);
	
//	for(int i=0; i<assignmentV.size(); i++){
//		cout << i << " " << assignmentV[i].chr << " "<< assignmentV[i].locus << " ";
//		for(int j=0; j<assignmentV[i].prob.size(); j++){
//			cout << assignmentV[i].prob[j]<< " ";
//		}
//		cout << endl;
//	}
	for(int i=0; i<assignmentV.size(); i++){
		for(int j=0; j<assignmentV.size(); j++){
			
			double sum=0.0;
			for(int k=0; k<assignmentV[i].prob.size(); k++){
				sum+=assignmentV[i].prob[k]*assignmentV[j].prob[k];
			}
			if(i==j){
				printf("%4.3f\t",1.0);
			}
			else{
				printf("%4.3f\t",sum);
			}
		}
		printf("\n");
	}
	
	return EXIT_SUCCESS;
	
}

void InputData(vector <assignment> &assignmentV, std::string infile1){
    
    ifstream file;
    file.open(infile1.c_str());
    string tmp;
	getline(file,tmp);
	vector <string> tmpV;
    trim(tmp);
//	cout << tmp << endl;
	getline(file,tmp);
//	cout << tmp << endl;
	
	getline(file,tmp);
	trim(tmp);
	boost::split(tmpV,tmp,boost::is_any_of("\t"));
//	cout << tmpV.size() << endl;
	file.close();
	int nClusters=tmpV.size()-2;
	
	file.open(infile1.c_str());
	getline(file,tmp);
	getline(file,tmp);
	
	while(!file.eof()){
		int i;
		assignment tmpA;
		file >> i;
		tmpA.chr=i;
		file >> i;
		tmpA.locus=i;
		double d;
		for(int k=0; k<nClusters; k++){
			file >>d;
			tmpA.prob.push_back(d);
		}
		if(!file.eof())
			assignmentV.push_back(tmpA);
	}
		
	file.close();
  // cout << "Read in data: Lines:" << assignmentV.size() << " Nclusters:" << assignmentV[0].prob.size() << endl;
    
}



