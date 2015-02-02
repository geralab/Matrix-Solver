/* Gerald Blake
 * CS 3513 Numerical Methods
 * 
 *
 *
 *
 */
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<cmath>
#include<vector>
using namespace std;
void openFile(ifstream &file,ofstream& outFile);
void readFile(ifstream &file, double** &Array, vector<vector<double>> &value,double* &b, double* &x, double *&temp, double* &rs,vector<double> &res,int &row, int &col, int &rowLength, int &colLength, char &type);
//void readFile(ifstream &file, double* &rowArray, int* &colIndex,double* &b,int &row, int &col, int &rowLength, int &colLength);
void printMatrix( double** Array,double *b,double* x, int rowLength, int colLength,ofstream &out);
void printMatrix(vector<vector<double>> &value, double *b,double *x, double *r,vector<double> res, int rowLength, int colLength, ofstream &out);
void NaiveGaussianElimination(double** &Array, double* &b,double* &x, int rowLength, int colLength, ofstream &out);
void GaussianEliminationPP(double** &Array, double* &b, double* &x, int rowLength, int colLength, ofstream &out);
void GuassSiedel(vector<vector<double>> &value, double* &b, double* &x,double* &rs,vector<double> &res, int rowLength, int colLength, double tol , double maxiter ,ofstream &out);
void Jacobi(vector<vector<double>> &value, double* &b, double* &x,double *&temp, double*& r,vector<double> &res, int rowLength, int colLength, double tol , double maxiter ,ofstream &out);
double residual(vector<vector<double>> &value, double* &b, double* &x,double * &rs,int rowLength, int colLength ,double tol,ofstream &out);
void SOR(vector<vector<double>> &value, double* &b, double* &x, double *&temp,double * &rs,vector<double> &res,int rowLength, int colLength, double tol , double maxiter, double &w,ofstream &out);
bool mArrayUsed = false;

int main()
{
	ifstream matrixFile;
	ofstream out;
	string fileName, matrixType;
	cout << setprecision(15);
	int row, col = 0;
	char type = ' ';
	int rowLength, colLength = 0;
	bool done = false;
	int task = 0;
	double **matrixArray = NULL;
	double *rowArray = NULL;
	vector<vector<double>> value;
	double *bVector = NULL;
	double *x = NULL;
	double *temp = NULL;
	double *r = NULL;
	double tol = .000000001;
	int maxiter = 10;
	double w;
	vector<double> res;
	// These statements refer to menu
	cout << "Press 0 to open a file and Populate Array" << endl;
	cout << "Press 1 to print Matrix " << endl;
	cout << "Press 2 to perform Guassian Elimination" << endl;
	cout << "Press 3 to perform Guassian Elimination w/ Partial Pivoting" << endl;
	cout << "Press 4 to perform the Guass Siedel Method"<< endl;
	cout << "Press 5 to perform the Jacobi Method" << endl;
	cout << "Press 6 to perform the SOR Method" << endl;
	cout << "Press 8 to repeat options" << endl;
	

	while(!done)
	{
		cout << "Which task would you like to perform?: ";
		cin >> task;
		switch(task)
		{
			case 0:
				openFile(matrixFile,out);
				readFile(matrixFile,matrixArray,value,bVector,x,temp,r,res,row,col,rowLength,colLength,type);
				break;
			case 1:
				if(type == 'd'){
					printMatrix(matrixArray,bVector,x,rowLength,colLength,out);
				}
				else if(type == 's')
				{
				    printMatrix(value,bVector,x,r,res,rowLength,colLength,out);
				}
				break;
			case 2:
				 NaiveGaussianElimination(matrixArray,bVector,x,rowLength,colLength,out);
				 
				break;
			case 3:
				 GaussianEliminationPP(matrixArray,bVector,x,rowLength,colLength,out);
				
				break;
			case 4:
				cout << "\nPlease enter the maximum amount of iterations and the tolerance: ";
				cin >> maxiter >> tol;
				GuassSiedel(value, bVector,x,r,res,rowLength,colLength,tol,maxiter,out);	
				break;
			case 5:
				cout << "\nPlease enter the maximum amount of iterations and the tolerance: ";
				cin >> maxiter >> tol;
				Jacobi(value, bVector,x,temp,r,res,rowLength,colLength,tol,maxiter,out);
				break;
				
			case 6:
				cout << "\nPlease enter the maximum amount of iterations, tolerance, and Omega(w): ";
				cin >> maxiter >> tol >> w;
				SOR(value,bVector,x,temp,r,res,rowLength,colLength,tol,maxiter,w,out);
				break;
			case 7:
				break;
			case 8:
				cout << "\n\nPress 0 to open a file and Populate Array" << endl;
				cout << "Press 1 to print Matrix " << endl;
				cout << "Press 2 to perform Guassian Elimination" << endl;
				cout << "Press 3 to perform Guassian Elimination w/ Partial Pivoting" << endl;
				cout << "Press 4 to perform the Guass Siedel Method"<< endl;
				cout << "Press 5 to perform the Jacobi Method" << endl;
				cout << "Press 6 to perform the SOR Method" << endl;
				cout << "Press 8 to repeat options" << endl;
				break;
			case 9:
				// deallocate memory
				if(mArrayUsed)
				{
					for(int i = 0; i < colLength; i++)
							delete [] matrixArray[i];
					delete [] matrixArray;
				}
				else 
				{
					delete [] matrixArray;
				}
				
				delete [] bVector;
				delete [] x;
				delete [] r;
				delete [] rowArray;
				delete [] temp;
				value.resize(0);
				value.clear();
				res.resize(0);
				res.clear();
				out.close();
				matrixFile.close();
				done = true;
				
		}
		cout << "\n\n";
	}

	return 0;
}

void openFile(ifstream &file, ofstream& outFile)
{
    string name;
	file.close();
	outFile.close();
	cout << "Enter file name: ";
	cin >> name;
	file.open(name.c_str());
}

void readFile(ifstream &file, double** &Array,vector<vector<double>> &value, double* &b, double* &x, double *&temp,double* &rs,vector<double> &res, int &row, int &col, int &rowLength, int &colLength, char &type)
{
	//get info on matrix
	file >> type >> rowLength >> colLength;
	//create dynamic array for rows
	
	type = tolower(type);
	if(type == 'd'){
	    mArrayUsed = true;
		Array = new double*[rowLength];
		b = new double[rowLength];
		x = new double[rowLength];
		//append columns to row positions
		for(int i = 0; i < colLength; i++)
			Array[i] = new double[colLength];

		//cout << rowLength << " " << colLength << endl;
		//populate array
		for (int row = 0; row < rowLength; row++)
		{
			for(int col = 0; col < colLength; col++)
			{
				file >> Array[row][col];
			}
		}

		for(int i = 0; i < rowLength; i++)
		{
			file >> b[i];
		}
	
	}
	else if(type == 's')
	{
		value.resize(0);
		res.resize(0);
		b = new double[rowLength];
		x = new double[rowLength];
		rs = new double[rowLength];
		temp = new double[rowLength];
		string line;
		getline(file,line);
		for(int r = 0; r < rowLength;r++){
			
			getline(file,line);
			int len = line.length() + 1;
			char *str = new char[len];
			strcpy(str,line.c_str());
			//cout << str << endl;
			char * pch;
			pch = strtok (str," ");
			value.push_back(vector<double>());
		    while (pch != NULL)
		    {
				double num = atof(pch);
				//cout << num <<endl;
				
				value[r].push_back(num);
				pch = strtok (NULL, " ");
		    }

			
		}
	
		for(int i = 0; i < rowLength; i++)
		{
			file >> b[i];
			x[i] = 0;
			temp[i] = 0;
			rs[i] = 0;
		}

		
	}
	
}

void printMatrix(double** Array, double *b, double* x, int rowLength, int colLength, ofstream &out)
{
	
	
	string outName;

	cout << "\n\nWhat would you like to print?" << endl;
	cout << "Press 0 to print matrix and solution " << endl;
	cout << "Press 1 to print matrix only" << endl;
	cout << "Press 2 to print solution only" << endl;
	cout << "Press 3 to repeat options" << endl;
	cout << "Press 4 to create a solution file" << endl;
	cout << "Press 5 to write to existing solution file" << endl;
	cout << "Press 6 to go back to main menu"<< endl;
	int option,option2 = 0;
	cin >> option;

	switch(option)
	{
		case 0:
			cout << " _____________________________SEPARATION BAR" << endl;
			for (int row = 0; row < rowLength; row++)
			{
				for(int col = 0; col < colLength; col++)
				{
					cout << Array[row][col] << " ";
				}
				//append b vector
				cout << "| " << b[row]<<"\n\n";
		
			}

				cout << "sol [";
				for(int i = 0; i < rowLength; i++)
				{
					cout << x[i] << ",";
				}
				cout << "]" << endl;
		break;

		case 1:
			cout << " _____________________________SEPARATION BAR" << endl;
			for (int row = 0; row < rowLength; row++)
			{
				for(int col = 0; col < colLength; col++)
				{
					cout << Array[row][col] << " ";
				}
				cout << "| " << b[row]<<"\n\n";
		
			}
			break;
		case 2:
			cout << " _____________________________SEPARATION BAR" << endl;
			cout << "sol [";
				for(int i = 0; i < rowLength; i++)
				{
					cout << x[i] << ",";
				}
				cout << "]" << endl;
				break;
		case 3:
			cout << "What would you like to print" << endl;
			cout << "Press 0 to print matrix and solution " << endl;
			cout << "Press 1 to print matrix only" << endl;
			cout << "Press 2 to print solution only" << endl;
			cout << "Press 3 to repeat options" << endl;
			cout << "Press 4 to create a solution file" << endl;
			cout << "Press 5 to write to existing file" << endl;
			cout << "Press 6 to go back to main menu"<< endl;
			cin >> option;
			
		case 4:
			cout << "\n\nWhat would you like to name your output file:";
			cin >> outName;
			cout << endl;
			out.open(outName.c_str());
			
		
			break;
		case 5:
			cout << "\n\nWhat would you like to print in your file?" << endl;
			cout << "Press 0 to print matrix and solution " << endl;
			cout << "Press 1 to print matrix only" << endl;
			cout << "Press 2 to print solution only" << endl;
			cout << "Press 3 to repeat options" << endl;
			cout << "Press 4 to go back to main menu"<< endl;
			cin >> option2;

			switch(option2)
			{
				case 0:
					out << " _____________________________SEPARATION BAR" << endl;
					for (int row = 0; row < rowLength; row++)
					{
						for(int col = 0; col < colLength; col++)
						{
							out << Array[row][col] << " ";
						}
						//append b vector
						out << "| " << b[row]<<"\n\n";
		
					}

						out << "sol [";
						for(int i = 0; i < rowLength; i++)
						{
							out << x[i] << ",";
						}
						out << "]" << endl;
				break;

				case 1:
					out << " _____________________________SEPARATION BAR" << endl;
					for (int row = 0; row < rowLength; row++)
					{
						for(int col = 0; col < colLength; col++)
						{
							out << Array[row][col] << " ";
						}
						out << "| " << b[row]<<"\n\n";
		
					}
					break;
				case 2:
					out << " _____________________________SEPARATION BAR" << endl;
					out << "sol [";
						for(int i = 0; i < rowLength; i++)
						{
							out << x[i] << ",";
						}
						out << "]" << endl;
						break;
				case 3:
					cout<< "File menu!!!" << endl;
					cout << "Press 0 to print matrix and solution " << endl;
					cout << "Press 1 to print matrix only" << endl;
					cout << "Press 2 to print solution only" << endl;
					cout << "Press 3 to repeat options" << endl;
					cout << "Press 4 to go back to main menu"<< endl;
					cin >> option2;
				case 4:
					break;
			}
				
				break;
		case 6:
			break;
	}
	
}
void printMatrix(vector<vector<double>> &value, double *b, double *x, double *r, vector<double> res, int rowLength, int colLength, ofstream &out)
{
	
	
	string outName;

	cout << "\n\nWhat would you like to print?" << endl;
	cout << "Press 0 to print matrix and solution " << endl;
	cout << "Press 1 to print matrix only" << endl;
	cout << "Press 2 to print solution only" << endl;
	cout << "Press 3 to repeat options" << endl;
	cout << "Press 4 to create a solution file" << endl;
	cout << "Press 5 to write to existing solution file" << endl;
	cout << "Press 6 to go back to main menu"<< endl;
	int option,option2 = 0;
	cin >> option;

	switch(option)
	{
		case 0:
			cout << " _____________________________SEPARATION BAR" << endl;
			for(unsigned int i = 0; i < value.size(); i++)
			{
				for(unsigned int j = 0; j < value[i].size();j++)
				{
					cout << value[i][j] << "  ";
				}
				
				cout << endl;
			}
			for(int i = 0 ; i < rowLength; i++)
			{
				cout << b[i] << endl;
			}

			cout << "sol [";
				for(int i = 0; i < rowLength; i++)
				{
					cout << x[i] << ",";
				}
				cout << "]\n" << endl;
				
				for(int i = 0; i < res.size(); i++)
				{
					cout <<i<< ":    " << res[i] << endl;;
				}
				

		break;

		case 1:
			cout << " _____________________________SEPARATION BAR" << endl;
			for(unsigned int i = 0; i < value.size(); i++)
			{
				for(unsigned int j = 0; j < value[i].size();j++)
				{
					cout << value[i][j] << "  ";
				}
				
				cout << endl;
			}
			for(int i = 0 ; i < rowLength; i++)
			{
				cout << b[i] << endl;
			}
			break;
		case 2:
			cout << "sol [";
				for(int i = 0; i < rowLength; i++)
				{
					cout << x[i] << ",";
				}
				cout << "]" << endl;
			
				break;
		case 3:
			cout << "What would you like to print" << endl;
			cout << "Press 0 to print matrix and solution " << endl;
			cout << "Press 1 to print matrix only" << endl;
			cout << "Press 2 to print solution only" << endl;
			cout << "Press 3 to repeat options" << endl;
			cout << "Press 4 to create a solution file" << endl;
			cout << "Press 5 to write to existing file" << endl;
			cout << "Press 6 to go back to main menu"<< endl;
			cin >> option;
			
		case 4:
			cout << "\n\nWhat would you like to name your output file:";
			cin >> outName;
			cout << endl;
			out.open(outName.c_str());
			
		
			break;
		case 5:
			cout << "\n\nWhat would you like to print in your file?" << endl;
			cout << "Press 0 to print matrix and solution " << endl;
			cout << "Press 1 to print matrix only" << endl;
			cout << "Press 2 to print solution only" << endl;
			cout << "Press 3 to repeat options" << endl;
			cout << "Press 4 to go back to main menu"<< endl;
			cin >> option2;

			switch(option2)
			{
				case 0:
					out << " _____________________________SEPARATION BAR" << endl;
				for(unsigned int i = 0; i < value.size(); i++)
				{
					for(unsigned int j = 0; j < value[i].size();j++)
					{
						out << value[i][j] << "  ";
					}
				
					out << endl;
				}
				for(int i = 0 ; i < rowLength; i++)
				{
					out << b[i] << endl;
				}

				out << "sol [";
					for(int i = 0; i < rowLength; i++)
					{
						out << x[i] << ",";
					}
					out << "]\n" << endl;
				
					for(int i = 0; i < res.size(); i++)
					{
						out <<i<< ":    " << res[i] << endl;;
					}
					
				break;

				case 1:
					out << " _____________________________SEPARATION BAR" << endl;
					for(unsigned int i = 0; i < value.size(); i++)
					{
						for(unsigned int j = 0; j < value[i].size();j++)
						{
							out << value[i][j] << "  ";
						}
				
						out << endl;
					}
					for(int i = 0 ; i < rowLength; i++)
					{
						out << b[i] << endl;
					}
					
					break;
				case 2:
					out << "sol [";
					for(int i = 0; i < rowLength; i++)
					{
						out << x[i] << ",";
					}
					out << "]\n" << endl;
				
					for(int i = 0; i < res.size(); i++)
					{
						out <<i<< ":    " << res[i] << endl;;
					}
						break;
				case 3:
					cout << "\n\nWhat would you like to print in your file?" << endl;
					cout << "Press 0 to print matrix and solution " << endl;
					cout << "Press 1 to print matrix only" << endl;
					cout << "Press 2 to print solution only" << endl;
					cout << "Press 3 to repeat options" << endl;
					cout << "Press 4 to go back to main menu"<< endl;
					cin >> option2;

					
				case 4:
					break;
			}
				
				break;
		case 6:
			break;
	}
	
}

void NaiveGaussianElimination(double** &Array, double* &b, double* &x, int rowLength, int colLength, ofstream &out)
{
	double mult, sum;
	bool pivotElementZero = false;
	for(int k = 0; k < rowLength - 1;k++)
	{
		for(int i = k + 1; i < colLength ; i++)
		{
			if( Array[i][k] != 0 )
			{
				if((Array[k][k]) == 0)
				{
					cout << "Pivot element is equal to zero, can't find solution using Naive Guassian method" << endl;
					if(out.is_open())
						out << "Pivot element is equal to zero, can't find solution using Naive Guassian method" << endl;

					return;
				}
				mult = (Array[i][k]) / (Array[k][k]);

				for(int j = k; j < colLength; j++)
				{
					Array[i][j] -= (mult * Array[k][j]);
					//printMatrix(Array,b,rowLength,colLength);
				}

				b[i] -= (mult * b[k]);
			}
		}
	}

	//backward substitution
		for(int i = rowLength - 1; i >= 0; i--)
		{
		   sum = 0.0;
		  for(int j = i + 1;j < colLength; j++)
		  {
		      sum += (Array[i][j] * x[j]);
		  }
		  if(Array[i][i] == 0)
		  {
			  cout << "Pivot element 0 and Rank Deficient - Naive Guassian" << endl;
			if(out.is_open())
			  out << "Pivot element 0 and Rank Deficient  - Naive Guassian" << endl;
			  return;
		  }
		   x[i] = (b[i] - sum) / (Array[i][i]);
		}

}

void GaussianEliminationPP(double** &Array, double* &b, double* &x, int rowLength, int colLength, ofstream &out)
{
	double mult, sum;
	double amax;
	int m;
	double temp;
	for(int k = 0; k < rowLength - 1;k++)
	{
		amax = abs(Array[k][k]);
		m = k;
		for (int j = k + 1; j < rowLength;j++)
		{
		   if(abs(Array[j][k]) > amax)
		   {
		      amax = abs(Array[j][k]);
		      m = j;
		   }
		}
		if( m != k) {
		   for(int j = k; j < rowLength; j++)
		   {
		      temp = Array[k][j];
		      Array[k][j] = Array[m][j];
		      Array[m][j] = temp;
		   }
		   temp = b[k];
		   b[k] = b[m];
		   b[m] = temp;
		}

		for(int i = k + 1; i < colLength ; i++)
		{
			if( Array[i][k] != 0 )
			{
				if((Array[k][k]) == 0)
				{
					cout << "Pivot element is equal to zero - Gaussian Partial Pivot" << endl;
					if(out.is_open())
						out << "Pivot element is equal to zero - Gaussian Partial Pivot" << endl;
					return;
				}
				mult = (Array[i][k]) / (Array[k][k]);

				for(int j = k; j < colLength; j++)
				{
					Array[i][j] -= (mult * Array[k][j]);
					
				}

				b[i] -= (mult * b[k]);
			}
		}
	}

	//backward substitution
		for(int i = rowLength - 1; i >= 0; i--)
		{
		   sum = 0.0;
		  for(int j = i + 1;j < colLength; j++)
		  {
		      sum += ((Array[i][j]) * (x[j]));
		  }

		  if(Array[i][i] == 0)
		  {
			  cout << "Pivot element 0 and Rank Deficient - Gaussian Partial Pivot" << endl;
			  if(out.is_open())
				out << "Pivot element is equal to zero - Gaussian Partial Pivot" << endl;
			  return;
		  }
		   x[i] = ((b[i]) - (sum)) / (Array[i][i]);
		}

}


void GuassSiedel(vector<vector<double>> &value, double* &b, double* &x, double* &rs,vector<double> &res, int rowLength, int colLength, double tol , double maxiter ,ofstream &out)
{
	
	double aij = 0;
	int j = 0;
	bool tolbool = false;
	for(int i = 0;(i < maxiter && !tolbool); i++)
	{
	
		for ( int r = 0; r < rowLength; r++)
		{
			double sum  = 0.0;
			double diag = 0;
			for(int c = 0; c  < value[r].size(); c+=2)
			{
				if(value[r][c] == r)
				{
					diag = value[r][c + 1];
				}
				else
				{
					j = (int)value[r][c];
					aij = value[r][c+1];
					sum += x[j] * aij;
				}
			}

			x[r] = (b[r] - sum) / diag;
	   }
		res.push_back(residual(value,b,x,rs,rowLength,colLength,tol ,out));
		if(res.size() > 0)
		{
			if(res[res.size() - 1] < tol)
				tolbool = true;

		}
		
	}

	
}

void Jacobi(vector<vector<double>> &value, double* &b, double* &x, double *&temp,double * &rs,vector<double> &res,int rowLength, int colLength, double tol , double maxiter ,ofstream &out)
{
	double aij;
	int j;
	bool tolbool = false;
	for(int i = 0;(i < maxiter && !tolbool); i++)
	{
		for ( int r = 0; r < rowLength; r++)
		{
			double sum  = 0.0;
			double diag = 0;
			for(int c = 0; c  < value[r].size(); c+=2)
			{
				if(value[r][c] == r)
				{
					diag = value[r][c + 1];
				}
				else
				{
					j = (int)value[r][c];
					aij = value[r][c+1];
					sum += x[j] * aij;
				}
			}

			temp[r] = (b[r] - sum) / diag;
	   }
		for(int k = 0; k < rowLength;k++)
		{
			x[k]  = temp[k];
		}
		res.push_back(residual(value,b,x,rs,rowLength,colLength,tol ,out));
		if(res.size() > 0)
		{
			if(res[res.size() - 1] < tol)
				tolbool = true;

		}
	}


}


double residual(vector<vector<double>> &value, double* &b, double* &x,double * &rs,int rowLength, int colLength , double tol, ofstream &out)
{
	double e = tol + 1;
	for(int r = 0 ; r < rowLength; r++)
	{
		double sum = 0.0;
		int j;
		double aij;
		for(int c = 0; c < value[r].size();c+=2)
		{

			j = (int)value[r][c];
			aij = value[r][c+1];
			sum += x[j] * aij;
			
		}
		rs[r] = b[r] - sum;
		e += rs[r]*rs[r];
	}
	e = sqrt(e);
	return e;
}


void SOR(vector<vector<double>> &value, double* &b, double* &x, double *&temp,double * &rs,vector<double> &res,int rowLength, int colLength, double tol , double maxiter, double &w ,ofstream &out)
{
	double aij = 0.0;
	int j = 0;
	double te;
	bool tolbool = false;
	for(int i = 0;(i < maxiter && !tolbool); i++)
	{
		for ( int r = 0; r < rowLength; r++)
		{
			double sum  = 0.0;
			double diag = 0;
			for(int c = 0; c  < value[r].size(); c+=2)
			{
				if(value[r][c] == r)
				{
					diag = value[r][c + 1];
				}
				else
				{
					j = (int)value[r][c];
					aij = value[r][c+1];
					sum += x[j] * aij;
				}
			}

			te = (b[r] - sum) / diag;
			x[r] = (1 - w) * x[r] + (w) * te;
	   }
		res.push_back(residual(value,b,x,rs,rowLength,colLength,tol ,out));
		if(res.size() > 0)
		{
			if(res[res.size() - 1] < tol)
				tolbool = true;

		}
		
	}


}