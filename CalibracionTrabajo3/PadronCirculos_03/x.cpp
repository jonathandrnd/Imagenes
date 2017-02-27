/* rename example */
#include <stdio.h>
#include <iostream>
#include <cstring>
#include <sstream>

using namespace std;
string tos(int t){
	stringstream st;
	st<<t;
	return st.str();
}

string f(int clase,int id){
	string dev="file_";
	dev+=tos(id);
	dev+=".jpg";
	return dev;
}

int main (){
	int result;
  	int clase=1;
  	int cont=0;
	for(int i=0;i<=1000;i++){
		string sold=f(clase,i);
		string snew=f(clase,cont);
  		result= rename( sold.c_str(), snew.c_str() );

  		if ( result == 0 ){
    		cout<<"yeeeeeee"<<endl;
			cont++;
		}else{
    		cout<<"no"<<endl;
		}
	}
  return 0;
}

