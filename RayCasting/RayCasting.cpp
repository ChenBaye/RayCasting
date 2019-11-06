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
#define SIZE	0.1	//单个像素的大小
#define INOUT(distance) (Radius-distance)>EPSILON?-1:(distance-Radius)>EPSILON?1:0	
//距球心distance的点是否在球内-1-内 1-外 0-球面

static float center[3] = { 0,0,0 };	//球心坐标

static float eye[3] = { -2,0,0 };	//视点坐标

static float image[3] = {-1,0,0};	//屏幕中心坐标，与x轴垂直

static float color_void[3] = { 0,0,0 };	//球体外的颜色

static float a_void = 0.005;	//球体外不透明度

static float color_sphere[3] = { 1,1,1 };	//球体的颜色

static float a_sphere = 0.015;	//球体不透明度




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

	float delta_x = position1[0] - position2[0];
	float delta_y = position1[1] - position2[1];
	float delta_z = position1[2] - position2[2];

	paraments[0] = 1;
	paraments[1] = delta_y / delta_x;
	paraments[2] = delta_z / delta_x;
}

//由上一个采样点start，计算出下一个采样点end
void GeneratePoint(float* start, float* end, float* paraments) {



}

//end是否超出包围盒范围
bool Inbox(float* end) {

}


//计算像素值
void GenerateRGBA(float* origin, float* RGBA, float* paraments) {
	float start[3] = { origin[0],origin[1],origin[2] };
	float end[3] = {0,0,0};

	float c_in[3] = { color_sphere[0], color_sphere[1] ,color_sphere[2] };
	float a_in = a_sphere;

	float c_now[3] = {0,0,0};
	float a_now = 0;

	// 依次计算出采样点
	while (true) {
		GeneratePoint(start, end, paraments);
		if (Inbox(end)) {
			Calculate(0)
		}
		else {
			break;
		}
	}


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
	float image_position[3] = { image[0],SIZE*WIDTH/2,SIZE*HEIGH/2};

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
	
	fclose(fp);
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