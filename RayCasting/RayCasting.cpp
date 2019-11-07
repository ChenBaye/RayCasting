#include <windows.h>  //Windows Header
#include <gl/gl.h>
#include <gl/glu.h>
#include <gl/glut.h>
#include <math.h>
#include <stdio.h>

#define EPSILON 0.000001	//浮点数比较大小的精度
#define cube_l	0.05   //cube的边长 
#define Radius	1		//单位球
#define LEN		(int)(Radius/cube_l + 1)		//包围盒半条边有多少个cube
#define NUM		2*LEN+1	//包围盒一条边有多少个点
#define MIN		(-(LEN * cube_l))	//包围盒最小坐标
#define MAX		+(LEN * cube_l)		//包围盒最大坐标
#define DIS		0.1		//采样距离
#define WIDTH	200	//屏幕宽度
#define HEIGH	200	//屏幕高度
#define SIZE	0.1		//单个像素的大小
#define INOUT(distance) (Radius-distance)>EPSILON?-1:(distance-Radius)>EPSILON?1:0	
//距球心distance的点是否在球内	-1-内 1-外 0-球面

static float center[3] = { 0,0,0 };	//球心坐标

static float eye[3] = { -2,0,0 };	//视点坐标

static float image[3] = {-1,0,0};	//屏幕中心坐标，与x轴垂直

static float color_void[3] = { 0,0,0 };	//球体外的颜色

static float a_void = 0.005;	//球体外不透明度

static float color_sphere[3] = { 1,1,1 };	//球体的颜色

static float a_sphere = 0.015;	//球体不透明度


//每次前进DIS，xyz的增量
float delta_x;
float delta_y;
float delta_z;


//体素结构体
typedef struct Voxel
{
	float x;
	float y;
	float z;
	float distance;
	int in_out;
}Voxel;

Voxel array[NUM][NUM][NUM];	//所有cube体素
float Image[WIDTH*HEIGH * 4];	//采用RGBA图像

//由体素某个方向的坐标，计算该体素在该方向上的索引
int Position2Index(float Position) {
	return (Position - MIN) / cube_l;
}


//单线性插值
float Linear_Interpolation(float value_0, float position_0, float value_1, float position_1, float position) {

	return value_0 + (value_1 - value_0) * (position - position_0) / (position_1 - position_0);

}

//是否需要插值
int NeedInterpolation(int(*Index)[3]) {
	int in = 0;
	int out = 0;
	int temp = 0;
	for (int i = 0; i < 8; i++) {
		temp = array[Index[i][0]][Index[i][1]][Index[i][2]].in_out;
		if (temp == -1) {

			in = 1;
		}
		else if (temp == 1) {
			out = 1;
		}
	}
	
	if (in != 0 && out != 0) {
		return 0;	//需要插值
	}
	else if (in != 0) {
		return -1;	//全部位于球内部(含球面)
	}	
	else {			
		return 1;	//全部位于球外部(含球面)
	}
}

//求两点间距离
float Distance(float a[3], float b[3])
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
		(a[1] - b[1]) * (a[1] - b[1]) +
		(a[2] - b[2]) * (a[2] - b[2]));
}

//距原点距离
float Distance_0(float x, float y, float z) {
	return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

//由两点，生成直线的参数方程
//	(x-x0)/m = (y-y0)/n = (z-z0)/p
//此处取m=1
void GenerateLine(float*position1, float*position2, float*paraments) {

	float sub_x = position2[0] - position1[0];
	float sub_y = position2[1] - position1[1];
	float sub_z = position2[2] - position1[2];

	paraments[0] = 1;
	paraments[1] = sub_y / sub_x;
	paraments[2] = sub_z / sub_x;

	//保存x,y,z增量
	float SUM_2 = paraments[0] * paraments[0] + paraments[1] * paraments[1] + paraments[2] * paraments[2];
	float temp = sqrt(SUM_2);

	delta_x = DIS * paraments[0] / temp;
	delta_y = DIS * paraments[1] / temp;
	delta_z = DIS * paraments[2] / temp;
}

//由上一个采样点start，计算出下一个采样点end
void GeneratePoint(float* start, float* end, float* paraments) {
	//	(x-x0)/m = (y-y0)/n = (z-z0)/p
	//	此处取m=1
	
	end[0] = start[0] + delta_x;
	end[1] = start[1] + delta_y;
	end[2] = start[2] + delta_z;
}

//end是否超出包围盒范围
bool Inbox(float* end) {
	for (int i = 0; i < sizeof(end)/sizeof(float); i++) {
		if (end[i]<MIN || end[i]>MAX) {
			return false;
		}
	}

	return true;
}

//位于包围盒与image之间
bool Beforebox(float* end) {
	if (end[0] > image[0] && end[0] < MIN) {
		return true;
	}
	else {
		return false;
	}
}

// 返回所在体元中沿某方向体素最小坐标
float GetMin(float x) {
	return ((int)(x / cube_l)-1)*cube_l;
}

// 返回所在体元中沿某方向体素最大坐标
float GetMax(float x) {
	return ((int)(x / cube_l))*cube_l;
}

//计算c_now，a_now
void Calculate_C_A(float*end, float*c_now, float*a_now) {
	//先求出8个体素坐标
	float x_min = GetMin(end[0]);
	float x_max = GetMax(end[0]);

	float y_min = GetMin(end[1]);
	float y_max = GetMax(end[1]);

	float z_min = GetMin(end[2]);
	float z_max = GetMax(end[2]);

	//存储八个体素的索引
	int Index[8][3] = {
		{ Position2Index(x_min), Position2Index(y_min), Position2Index(z_min)},
		{ Position2Index(x_max), Position2Index(y_min), Position2Index(z_min)},
		{ Position2Index(x_max), Position2Index(y_max), Position2Index(z_min)},
		{ Position2Index(x_min), Position2Index(y_max), Position2Index(z_min)},
		{ Position2Index(x_min), Position2Index(y_min), Position2Index(z_max)},
		{ Position2Index(x_max), Position2Index(y_min), Position2Index(z_max)},
		{ Position2Index(x_max), Position2Index(y_max), Position2Index(z_max)},
		{ Position2Index(x_min), Position2Index(y_max), Position2Index(z_max)},
	};
	
	//是否需要插值
	int result = NeedInterpolation(Index);
	switch (result) {
		case -1:		//全在内部
			c_now[0] = color_sphere[0];
			c_now[1] = color_sphere[1];
			c_now[2] = color_sphere[2];
			*a_now = a_sphere;
			break;
		case 1:		//全在外部
			c_now[0] = color_void[0];
			c_now[1] = color_void[1];
			c_now[2] = color_void[2];
			*a_now = a_void;
			break;
		case 0:		//有内有外，需要插值
			//将三线性插值分为七次单线性插值
			//先取过end点垂直于x轴的平面截取体元得四个点
			float c_p1[3], a_p1;
			float c_p2[3], a_p2;
			float c_p3[3], a_p3;
			float c_p4[3], a_p4;

			for (int i = 0; i < 3; i++) {
				c_p1[0] = Linear_Interpolation();
			}

			float c_p5[3], a_p5;
			float c_p6[3], a_p6;
			float c_p7[3], a_p7;



			break;
	}

}


//计算像素值
void GenerateRGBA(float* origin, float* RGBA, float* paraments) {
	float start[3] = { origin[0],origin[1],origin[2] };
	float end[3] = {0,0,0};

	//经过采样点前
	float c_in[3] = { 0, 0 ,0 };
	float a_in = 0;

	//采样点
	float c_now[3] = {0,0,0};
	float a_now = 0;

	//经过采样点后
	float c_out[3] = { 0,0,0 };
	float a_out = 0;

	while (true) {
		// 计算出采样点
		GeneratePoint(start, end, paraments);

		if ((Beforebox(end) || Inbox(end)) && a_out<1) {		//不能超出包围盒，且累计不透明度不能超过1，但可以位于包围盒与image之间
			
			Calculate_C_A(end, c_now, &a_now);	//计算c_now，a_now

			a_out = a_in + a_now*(1-a_in);	//不透明度A
			c_out[0] = c_in[0]*a_in + c_now[0]*a_now*(1 - a_in);	//颜色R
			c_out[1] = c_in[1]*a_in + c_now[1]*a_now*(1 - a_in);	//颜色G
			c_out[2] = c_in[2]*a_in + c_now[2]*a_now*(1 - a_in);	//颜色B

			a_in = a_out;
			c_in[0] = c_out[0];
			c_in[1] = c_out[1];
			c_in[2] = c_out[2];
		}
		else {
			break;
		}
	}

	RGBA[0] = c_out[0];
	RGBA[1] = c_out[1];
	RGBA[2] = c_out[2];
	RGBA[3] = a_out;
}



//遍历包围盒，计算体素
void calculate_voxel() {
	float x = MIN;
	float y = MIN;
	float z = MIN;
	//计算体素
	for (int i = 0; i < NUM; i++, x += cube_l) {	//循环NUM次
		y = MIN;
		for (int j = 0; j < NUM; j++, y += cube_l) {
			z = MIN;
			for (int k = 0; k < NUM; k++, z += cube_l) {
				array[i][j][k].x = x;	//坐标
				array[i][j][k].y = y;
				array[i][j][k].z = z;
				array[i][j][k].distance = Distance_0(x, y, z);	//距球心距离
				array[i][j][k].in_out = INOUT(array[i][j][k].distance);	//球内还是球外
			}
		}
	}
}

//从前到后等距采样,计算颜色和不透明度
void RayCasting() {

	int count = 0;

	//像素取值
	float RGBA[4] = { 0,0,0,0 };

	//直线参数
	float paraments[3] = { 1,0,0 };

	//像素坐标
	float image_position[3] = { image[0],-SIZE*WIDTH/2,-SIZE*HEIGH/2};

	//遍历像素点
	for (int i = 0; i < HEIGH; i++) {
		for (int j = 0; j < WIDTH; j++) {
			//计算直线参数方程
			GenerateLine(eye, image_position, paraments);

			//计算像素值
			GenerateRGBA(image_position, RGBA, paraments);

			//拷贝像素
			Image[count+0] = RGBA[0];
			Image[count+1] = RGBA[1];
			Image[count+2] = RGBA[2];
			Image[count+3] = RGBA[3];
			count += 4;

			//下一像素坐标
			image_position[1] += SIZE;
			image_position[2] += SIZE;
		}
	}
}




//使用opengl展示
void display_voxel() {

	glClear(GL_COLOR_BUFFER_BIT);//清颜色缓冲区
	glVertex2f(1, 1);
	glBegin(GL_POINTS);
	

			//glVertex3f(x, y, z);
	
	glEnd();
	glFlush();
}


int main(int argc, char** argv) {

	calculate_voxel();






	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(800, 800);
	glutCreateWindow("单位球");
	glutDisplayFunc(display_voxel);
	glutMainLoop();
	int i = 0;
	scanf("%d", &i);
	scanf("%d", &i);

	return 0;
}