#include<iostream>
#include <opencv2/highgui.hpp>
#include <opencv/cv.h>
#include <fstream>
#include <Eigen/Eigen>
using namespace std;
using namespace cv;
using namespace Eigen;
class tdpants
{
public:
	double width[5];
	double height[2];
	int resolutionx;
	int resolutiony[2];
	int numberoffaces;
	Matrix<double, 60, 250> matx;
	Matrix<double, 60, 250> maty;
	Vector2d point2d[7];
	tdpants(double* p1, double* p2, int res)
	{
		for (int i = 0; i < 5; i++)
		{
			width[i] = p1[i];
		}
		height[0] = p2[0];
		height[1] = p2[1];
		resolutionx = res;
	}
	void printfortest()
	{
		cout << "width(0~4):";
		for (int i = 0; i < 5; i++)
		{
			cout<<width[i]<<" ";
		}
		cout << endl;
		cout << "height(0~1):";
		cout << height[0] << " " << height[1] << endl;
		cout << "resolution:" << resolutionx << endl;
	}
	void computepoint()
	{
		point2d[0] << 0, 0;
		point2d[1] << width[0], 0;
		point2d[2] << 0, height[0];
		point2d[3] << width[1], height[0];
		point2d[4] << width[2], height[0];
		point2d[5] << width[3], height[1];
		point2d[6] << width[4], height[1];
	}
	void printpoint()
	{
		for (int i = 0; i < 7; i++)
		{
			cout << "(" << point2d[i][0] << "," << point2d[i][1] << ")" << endl;
		}
	}
	void makematrix()
	{
		double dx = width[0] / resolutionx;
		double dy = dx;
		resolutiony[0] = height[0] / dy;
		resolutiony[1] = (height[1]-height[0]) / dy;
		//cout << resolutiony[0] << " " << resolutiony[1] << endl;
		double xinit[1000];
		double dxinit[1000];
		double k1 = (width[2] - width[0])*1.0 / height[0];
		double k2 = (width[3] - width[1])*1.0 / (height[1] - height[0]);
		double k3 = (width[4] - width[2])*1.0 / (height[1] - height[0]);
		//cout << k1 << " " << k2 << " " << k3 << endl;
		for (int i = 0; i <= resolutiony[0] + resolutiony[1]; i++)
		{
			if (i <= resolutiony[0])
			{
				xinit[i] = 0;
				dxinit[i] = (width[0] + i * dy * k1) / resolutionx;
			}
			else
			{
				int modifiedi = i - resolutiony[0];
				xinit[i] = width[1] + modifiedi * dy * k2;
				dxinit[i] = (width[2] - width[1] + modifiedi * dy * k3- modifiedi * dy * k2)/resolutionx;
			}
			//cout << xinit[i] << " " << dxinit[i] << endl;
		}
		for (int j = 0; j <= resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = 0; i <= resolutionx; i++)
			{
				matx(i, j) = xinit[j]+i*dxinit[j];
				maty(i, j) = j * dy;
			}
		}
		matx = -matx;
		maty = -maty;
	}
	void makexyz(double x = 0, double y = 0,double z=0)
	{
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\pant.xyz");
		for (int j = 0; j <= resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = resolutionx; i >= 0; i--)
			{
				ofile << -matx(i, j)+x << " " << maty(i, j)+y << " " << z << endl;
			}
			for (int i = 0; i <= resolutionx; i++)
			{
				ofile << matx(i, j)+x << " " << maty(i, j)+y << " " << z << endl;
			}
		}
		ofile.close();
	}
	int makemesh(double x = 0, double y = 0, double z = 0)
	{
		int indexf=0;
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\pant1.obj");
		for (int j = 0; j <= resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = resolutionx; i >= 0; i--)
			{
				ofile << "v " << -matx(i, j)+x << " " << maty(i, j)+y << " " << z << endl;
			}
			for (int i = 0; i <= resolutionx; i++)
			{
				ofile << "v " << matx(i, j)+x << " " << maty(i, j)+y << " " << z << endl;
			}
		}
		int numberperline = 2 * resolutionx + 2;
		for (int j = 0; j < resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = 0; i < numberperline-1; i++)
			{
				if (i == resolutionx)
				{
					continue;
				}
				ofile << "f " << j * numberperline + i+1 << " ";
				ofile << j * numberperline + i + 1+1 << " ";
				ofile << (j + 1)*numberperline + i + 1+1 << " ";
				ofile << (j + 1)*numberperline + i+1  << endl;
				indexf++;
			}
		}
		numberoffaces = indexf;
		ofile.close();
		return indexf;
	}

	int makeanothermesh(double x = 0, double y = 0, double z = 0)
	{
		int indexf = 0;
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\pant2.obj");
		for (int j = 0; j <= resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = resolutionx; i >= 0; i--)
			{
				ofile << "v " << matx(i, j) - x << " " << maty(i, j) + y << " " << z << endl;
			}
			for (int i = 0; i <= resolutionx; i++)
			{
				ofile << "v " << -matx(i, j) - x << " " << maty(i, j) + y << " " << z << endl;
			}
		}
		int numberperline = 2 * resolutionx + 2;
		for (int j = 0; j < resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = 0; i < numberperline - 1; i++)
			{
				if (i == resolutionx)
				{
					continue;
				}
				ofile << "f " << j * numberperline + i + 1 << " ";
				ofile << j * numberperline + i + 1 + 1 << " ";
				ofile << (j + 1)*numberperline + i + 1 + 1 << " ";
				ofile << (j + 1)*numberperline + i + 1 << endl;
				indexf++;
			}
		}
		numberoffaces = indexf;
		ofile.close();
		return indexf;
	}

	int makemeshdata(double v[][3], int f[][4],int pointoffset=0,double x = 0, double y = 0, double z = 0)
	{
		int indexv = 0;
		int indexf = 0;
		for (int j = 0; j <= resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = resolutionx; i >= 0; i--)
			{
				v[indexv][0] = -matx(i, j) + x;
				v[indexv][1] = maty(i, j) + y;
				v[indexv][2] = z;
				indexv++;
				
			}
			for (int i = 0; i <= resolutionx; i++)
			{
				v[indexv][0] = matx(i, j) + x;
				v[indexv][1] = maty(i, j) + y;
				v[indexv][2] = z;
				indexv++;
			}
		}
		int numberperline = 2 * resolutionx + 2;
		for (int j = 0; j < resolutiony[0] + resolutiony[1]; j++)
		{
			for (int i = 0; i < numberperline - 1; i++)
			{
				if (i == resolutionx)
				{
					continue;
				}
				f[indexf][0] = j * numberperline + i + 1 + pointoffset;
				f[indexf][1] = j * numberperline + i + 1 + 1 + pointoffset;
				f[indexf][2] = (j + 1)*numberperline + i + 1 + 1 + pointoffset;
				f[indexf][3] = (j + 1)*numberperline + i + 1 + pointoffset;
				indexf++;
			}
		}
		numberoffaces = indexf;
		return indexf;
	}
};

void readsmpl(double v[6890][3],int f[13776][3])
{
	string data;
	ifstream ifile;
	ifile.open("F:\\ClothSimulation\\data\\smpl0.obj");
	assert(ifile.is_open());
	for (int i = 0; i < 6890; i++)
	{
		//cout << i << endl;
		string data;
		getline(ifile, data);
		if (data[0] != 'v' || data[1] != ' ')
		{
			i--;
			continue;
		}
		int j = 2;
		char num1[3][20];
		int numflag[3] = { 0,0,0 };
		for (int k = 0; k < 3; k++)
		{
			while (data[j] != '\0'&&data[j] != ' '&&data[j] != '\n')
			{
				num1[k][numflag[k]] = data[j];
				numflag[k]++;
				j++;
			}
			j++;
			num1[k][numflag[k]] = '\0';
			stringstream temp;
			temp << num1[k];
			temp >> v[i][k];
			//cout << v[i][k] << " ";
		}
		//cout << endl;
	}
	for (int i = 0; i < 13776; i++)
	{
		//cout << i << endl;
		string data;
		getline(ifile, data);
		if (data[0] != 'f' || data[1] != ' ')
		{
			i--;
			continue;
		}
		int j = 2;
		char num1[3][20];
		int numflag[3] = { 0,0,0 };
		for (int k = 0; k < 3; k++)
		{
			while (data[j] != '\0'&&data[j] != ' '&&data[j] != '\n')
			{
				num1[k][numflag[k]] = data[j];
				numflag[k]++;
				j++;
			}
			j++;
			num1[k][numflag[k]] = '\0';
			stringstream temp;
			temp << num1[k];
			temp >> f[i][k];
			//cout << f[i][k] << " ";
		}
		//cout << endl;
		//system("pause");
	}
	ifile.close();
}

void toghthertomesh(double pant_v[1188][3],int pant_f[2120][4], double smpl_v[6890][3],int smpl_f[13776][3], int pantfaces,int pointoffset,double x=0,double y=0,double z=0)
{
	ofstream ofile;
	ofile.open("C:\\Users\\lenovo\\Desktop\\data\\pantall.obj");
	for (int i = 0; i < 6890; i++)
	{
		ofile << "v " << smpl_v[i][0] << " " << smpl_v[i][1] << " " << smpl_v[i][2] << endl;
	}
	for (int i = 0; i < 1188; i++)
	{
		ofile << "v " << pant_v[i][0] << " " << pant_v[i][1] << " " << pant_v[i][2] << endl;
	}
	for (int i = 0; i < 1188; i++)
	{
		ofile << "v " << -pant_v[i][0]-x << " " << pant_v[i][1]+y << " " << pant_v[i][2]+z << endl;
	}
	for (int i = 0; i < 13776; i++)
	{
		ofile << "f " << smpl_f[i][0] << " " << smpl_f[i][1] << " " << smpl_f[i][2] << endl;
	}
	for (int i = 0; i < pantfaces; i++)
	{
		ofile << "f " << pant_f[i][0] << " " << pant_f[i][1] << " " << pant_f[i][2] << " " << pant_f[i][3] << endl;
	}
	for (int i = 0; i < pantfaces; i++)
	{
		ofile << "f " << pant_f[i][0]+pointoffset << " " << pant_f[i][1] + pointoffset << " " << pant_f[i][2] + pointoffset << " " << pant_f[i][3] + pointoffset << endl;
	}
	ofile.close();
}

void modifyparameter(double p1[5], double p2[2],double size)
{
	for (int i = 0; i < 5; i++)
	{
		p1[i] /= size;
	}
	for (int i = 0; i < 2; i++)
	{
		p2[i] /= size;
	}
}


int main()
{
	double smpl_v[6890][3];
	int smpl_f[13776][3];
	readsmpl(smpl_v, smpl_f);

	double p1[5] = { 7,-4,7.5,0.6,6.1 };
	double p2[2] = { 10,21 };
	modifyparameter(p1, p2, 23);
	int res = 20;
	tdpants PANT = tdpants(p1,p2,res);
	//PANT.printfortest();
	PANT.computepoint();
	//PANT.printpoint();
	PANT.makematrix();
	//PANT.makexyz(0,0,0);
	PANT.makemesh(0, -0.17, 0.2);
	double pant_v[1188][3];
	int pant_f[2120][4];
	PANT.makemeshdata(pant_v,pant_f,6890,0,-0.17,0.18);

	//toghthertomesh(pant_v, pant_f, smpl_v, smpl_f,PANT.numberoffaces,1188,0,0,-0.32);
	PANT.makeanothermesh(0, -0.17, -0.14);
	system("D:\\MeshLab\\meshlab.exe F:\\ClothSimulation\\data\\pant1.obj");
	return 0;
}