#include<iostream>
#include <fstream>
#include<vector>
#include <assert.h>
#include<string>
#include <sstream>
using namespace std;

class point2d
{
public:
	double x;
	double y;
};

class tdshirt
{
public:
	double width[4];
	double height[3];
	int resolutionx[3];
	int resolutiony[3];
	vector<vector<point2d>> point;
	int pointnumber;
	int facenumber;
	tdshirt(double* w, double* h, int rx0)
	{
		width[0] = w[0];
		width[1] = w[1];
		width[2] = w[2];
		width[3] = w[3];
		height[0] = h[0];
		height[1] = h[1];
		height[2] = h[2];
		resolutionx[0] = rx0;
		resolutionx[1] = int(rx0*width[1] / width[0]);
		resolutionx[2] = int(rx0*width[2] / width[0]);
		resolutiony[0] = int(rx0*height[0] / width[0]);
		resolutiony[1] = int(rx0*height[1] / width[0]);
		resolutiony[2] = int(rx0*height[2] / width[0]);
		point.resize(resolutiony[0] + resolutiony[1] + resolutiony[2] + 1);
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			point[y].reserve(1000);

		}
	}
	void printparameter()
	{
		cout << width[0] << " " << width[1] << " " << width[2] << " " << width[3] << endl;
		cout << height[0] << " " << height[1] << " " << height[2] << endl;
		cout << resolutionx[0] << " " << resolutionx[1] << " " << resolutionx[2] << endl;
		cout << resolutiony[0] << " " << resolutiony[1] << " " << resolutiony[2] << endl;
	}
	int makepoint()
	{
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			if (y <= resolutiony[0])
			{
				point[y].resize(resolutionx[0] + 1);
				for (int x = 0; x <= resolutionx[0]; x++)
				{
					point[y][x].x = width[0] / resolutionx[0] * x;
					point[y][x].y = height[0] / resolutiony[0] * y;
				}
			}
			else
				if (y <= resolutiony[0] + resolutiony[1])
				{
					point[y].resize(resolutionx[0] + resolutionx[1] + resolutionx[2] + 1);
					for (int x = 0; x <= resolutionx[0]; x++)
					{
						point[y][x].x = width[0] / resolutionx[0] * x;
						point[y][x].y = height[0] + height[1] / resolutiony[1] * (y - resolutiony[0]);
					}
				}
				else
				{
					point[y].resize(resolutionx[0] + resolutionx[1] + 1);
					for (int x = 0; x <= resolutionx[0]; x++)
					{
						point[y][x].x = width[0] / resolutionx[0] * x;
						point[y][x].y = height[0] + height[1] + height[2] / resolutiony[2] * (y - resolutiony[0] - resolutiony[1]);
					}
				}
		}
		for (int y = resolutiony[0]+1; y <= resolutiony[0] + resolutiony[1]; y++)
		{
			for (int x = resolutionx[0] + 1; x <= resolutionx[0] + resolutionx[1]; x++)
			{
				point[y][x].x = width[0] + width[1] / resolutionx[1] * (x - resolutionx[0]);
				point[y][x].y = height[0] + height[1] / resolutiony[1] * (y - resolutiony[0]);
			}
			for (int x = resolutionx[0] + resolutionx[1] + 1; x <= resolutionx[0] + resolutionx[1] + resolutionx[2]; x++)
			{
				point[y][x].x = width[0] + width[1] + width[2] / resolutionx[2] * (x - resolutionx[0] - resolutionx[1]);
				point[y][x].y = height[0] + height[1] / resolutiony[1] * (y - resolutiony[0]);
			}
		}
		double k = (width[3] - width[1]) / height[2];
		for (int y = resolutiony[0] + resolutiony[1]+1; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			for (int x = resolutionx[0] + 1; x <= resolutionx[0] + resolutionx[1]; x++)
			{
				point[y][x].x = width[0] + (width[1] + k * (y - resolutiony[0] - resolutiony[1])) / resolutionx[1] * (x - resolutionx[0]);
				point[y][x].y = height[0] + height[1] + height[2] / resolutiony[2] * (y - resolutiony[0] - resolutiony[1]);
			}
		}
		return 0;
	}
	void shiftxy()
	{
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = 0;
			while (x < point[y].size())
			{
				point[y][x].x *= -1;
				point[y][x].y *= -1;
				x++;
			}
		}
	}
	void makexyz(double x_offset = 0, double y_offset = 0, double z_offset = 0)
	{
		pointnumber = 0;
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\shirt.xyz");
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = point[y].size() - 1;
			while (x >= 0)
			{
				ofile << "v " << -point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x--;
				pointnumber++;
			}
			x = 0;
			while (x < point[y].size())
			{
				ofile << "v " <<point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x++;
				pointnumber++;
			}
		}
		ofile.close();
	}
	void makemesh(int pointoffset = 0, double x_offset = 0, double y_offset = 0, double z_offset = 0)
	{
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\shirt1.obj");
		pointnumber = 0;
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = point[y].size()-1;
			while (x>=0)
			{
				ofile << "v "<<-point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x--;
				pointnumber++;
			}
			x = 0;
			while (x < point[y].size())
			{
				ofile << "v " << point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x++;
				pointnumber++;
			}
		}
		int numberperline =0;
		int totalnumber = pointoffset;
		facenumber = 0;
		numberperline = 2 * resolutionx[0] + 2;
		for (int y = 0; y < resolutiony[0] ; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber+numberperline + x+1 << " ";
				ofile << totalnumber +numberperline + x << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[1]+resolutionx[2]+x+1 << " ";
			ofile << totalnumber + numberperline + resolutionx[1]+resolutionx[2]+x << endl;
			facenumber++;
		}
		totalnumber += numberperline;
		numberperline = 2 * (resolutionx[0] + resolutionx[1] + resolutionx[2]) + 2;
		for (int y = 0; y < resolutiony[1]-1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber +x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber + numberperline + x+1 << " ";
				ofile << totalnumber + numberperline + x  << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		totalnumber += resolutionx[2];
		numberperline = 2 * (resolutionx[0] + resolutionx[1]) + 2;
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + x + 1 << " ";
			ofile << totalnumber + numberperline+resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x << endl;
			facenumber++;
		}
		totalnumber += 2*resolutionx[0]+2*resolutionx[1]+resolutionx[2]+2;
		for (int y = 0; y < resolutiony[2]-1; y++) 
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber +x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber + numberperline + x+1 << " ";
				ofile << totalnumber + numberperline + x  << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		ofile.close();
	}
	void makeanothermesh(int pointoffset=0,double x_offset = 0, double y_offset = 0, double z_offset = 0)
	{
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\shirt2.obj");
		pointnumber = 0;
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = point[y].size() - 1;
			while (x >= 0)
			{
				ofile << "v " << -point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x--;
				pointnumber++;
			}
			x = 0;
			while (x < point[y].size())
			{
				ofile << "v " << point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset << endl;
				x++;
				pointnumber++;
			}
		}
		int numberperline = 0;
		int totalnumber = pointoffset;
		facenumber = 0;
		numberperline = 2 * resolutionx[0] + 2;
		for (int y = 0; y < resolutiony[0]; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x <<"  " ;
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + x + 1 << endl;
			facenumber++;
		}
		totalnumber += numberperline;
		numberperline = 2 * (resolutionx[0] + resolutionx[1] + resolutionx[2]) + 2;
		for (int y = 0; y < resolutiony[1] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		totalnumber += resolutionx[2];
		numberperline = 2 * (resolutionx[0] + resolutionx[1]) + 2;
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + x + 1 << endl;
			facenumber++;
		}
		totalnumber += 2 * resolutionx[0] + 2 * resolutionx[1] + resolutionx[2] + 2;
		for (int y = 0; y < resolutiony[2] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
				facenumber++;
			}
			totalnumber += numberperline;
		}
		ofile.close();
	}
	void meshwithsmpl(double v[6890][3], int f[13776][3],double x_offset=0,double y_offset=0,double z_offset1=0.2,double z_offset2=-0.2)
	{
		ofstream ofile;
		ofile.open("F:\\ClothSimulation\\data\\shirtall.obj");
		for (int i = 0; i < 6890; i++)
		{
			ofile << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << endl;
		}
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = point[y].size() - 1;
			while (x >= 0)
			{
				ofile << "v " << -point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset1 << endl;
				x--;
			}
			x = 0;
			while (x < point[y].size())
			{
				ofile << "v "<<point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset1 << endl;
				x++;
			}
		}
		for (int y = 0; y <= resolutiony[0] + resolutiony[1] + resolutiony[2]; y++)
		{
			int x = point[y].size() - 1;
			while (x >= 0)
			{
				ofile << "v " << -point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset2 << endl;
				x--;
			}
			x = 0;
			while (x < point[y].size())
			{
				ofile << "v "<<point[y][x].x + x_offset << " " << point[y][x].y + y_offset << " " << z_offset2 << endl;
				x++;
			}
		}
		for (int i = 0; i < 13776; i++)
		{
			ofile << "f " << f[i][0] << " " << f[i][1] << " " << f[i][2] << endl;
		}
		int numberperline = 0;
		int totalnumber = 6890;
		numberperline = 2 * resolutionx[0] + 2;
		for (int y = 0; y < resolutiony[0]; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + numberperline + x << endl;
			}
			totalnumber += numberperline;
		}
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x << endl;
		}
		totalnumber += numberperline;
		numberperline = 2 * (resolutionx[0] + resolutionx[1] + resolutionx[2]) + 2;
		for (int y = 0; y < resolutiony[1] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + numberperline + x << endl;
			}
			totalnumber += numberperline;
		}
		totalnumber += resolutionx[2];
		numberperline = 2 * (resolutionx[0] + resolutionx[1]) + 2;
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x << endl;
		}
		totalnumber += 2 * resolutionx[0] + 2 * resolutionx[1] + resolutionx[2] + 2;
		for (int y = 0; y < resolutiony[2] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + x + 1 << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + numberperline + x << endl;
			}
			totalnumber += numberperline;
		}
		numberperline = 0;
		totalnumber = 6890+pointnumber;
		numberperline = 2 * resolutionx[0] + 2;
		for (int y = 0; y < resolutiony[0]; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x << "  ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
			}
			totalnumber += numberperline;
		}
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x << " ";
			ofile << totalnumber + numberperline + resolutionx[1] + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + x + 1 << endl;
		}
		totalnumber += numberperline;
		numberperline = 2 * (resolutionx[0] + resolutionx[1] + resolutionx[2]) + 2;
		for (int y = 0; y < resolutiony[1] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
			}
			totalnumber += numberperline;
		}
		totalnumber += resolutionx[2];
		numberperline = 2 * (resolutionx[0] + resolutionx[1]) + 2;
		for (int x = 1; x < numberperline; x++)
		{
			ofile << "f " << totalnumber + x << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x << " ";
			ofile << totalnumber + numberperline + resolutionx[2] + x + 1 << " ";
			ofile << totalnumber + x + 1 << endl;
		}
		totalnumber += 2 * resolutionx[0] + 2 * resolutionx[1] + resolutionx[2] + 2;
		for (int y = 0; y < resolutiony[2] - 1; y++)
		{
			for (int x = 1; x < numberperline; x++)
			{
				ofile << "f " << totalnumber + x << " ";
				ofile << totalnumber + numberperline + x << " ";
				ofile << totalnumber + numberperline + x + 1 << " ";
				ofile << totalnumber + x + 1 << endl;
			}
			totalnumber += numberperline;
		}
		ofile.close();
	}
};

void readsmpl(double v[6890][3], int f[13776][3])
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

void modifyparameter(double p1[4], double p2[3], double size)
{
	for (int i = 0; i < 4; i++)
	{
		p1[i] /= size;
	}
	for (int i = 0; i < 3; i++)
	{
		p2[i] /= size;
	}
}

int main()
{
	double width[4] = { 1.0, 3.5, 7.0,3.51 };
	double height[3] = { 0.3, 3.0, 8.0 };
	modifyparameter(width, height, 17.0);
	double x_offset = 0;
	double y_offset = 0.34;
	double z_offset1 = 0.3;
	double z_offset2 = -0.3;
	int rx0 = 5;
	tdshirt SHIRT = tdshirt(width, height, rx0);
	SHIRT.printparameter();
	SHIRT.makepoint();
	SHIRT.shiftxy();
	SHIRT.makemesh(0,x_offset,y_offset,z_offset1);
	SHIRT.makeanothermesh(0, x_offset, y_offset, z_offset2);
	double v[6890][3];
	int f[13776][3];
	readsmpl(v, f);
	SHIRT.meshwithsmpl(v,f, x_offset, y_offset, z_offset1, z_offset2);
	//cout<<"pointnumber="<<SHIRT.pointnumber<<endl;
	//cout << "facenumber=" << SHIRT.facenumber<<endl;
	system("D:\\MeshLab\\meshlab.exe F:\\ClothSimulation\\data\\shirtall.obj");
	return 0;
}