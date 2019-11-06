#include <windows.h>  //Windows Header
#include <gl/gl.h>
#include <gl/glu.h>
#include <gl/glut.h>
#include <math.h>
#include <stdio.h>

#define EPSILON 0.000001	//�������Ƚϴ�С�ľ���
#define cube_l	0.05   //cube�ı߳� 
#define Radius	1		//��λ��
#define LEN		(int)(Radius/cube_l + 1)		//��Χ�а������ж��ٸ�cube
#define NUM		2*LEN+1	//��Χ��һ�����ж��ٸ���
#define MIN		(-(LEN * cube_l))	//��Χ����С����
#define MAX		+(LEN * cube_l)		//��Χ���������
#define DIS		0.1		//��������
#define WIDTH	200	//��Ļ���
#define HEIGH	200	//��Ļ�߶�
#define SIZE	0.1	//�������صĴ�С
#define INOUT(distance) (Radius-distance)>EPSILON?-1:(distance-Radius)>EPSILON?1:0	
//������distance�ĵ��Ƿ�������-1-�� 1-�� 0-����

static float center[3] = { 0,0,0 };	//��������

static float eye[3] = { -2,0,0 };	//�ӵ�����

static float image[3] = {-1,0,0};	//��Ļ�������꣬��x�ᴹֱ

static float color_void[3] = { 0,0,0 };	//���������ɫ

static float a_void = 0.005;	//�����ⲻ͸����

static float color_sphere[3] = { 1,1,1 };	//�������ɫ

static float a_sphere = 0.015;	//���岻͸����




//���ؽṹ��
typedef struct Voxel
{
	float x;
	float y;
	float z;
	float distance;
	int in_out;
}Voxel;

Voxel array[NUM][NUM][NUM];	//����cube����
float Image[WIDTH*HEIGH * 4];	//����RGBAͼ��

//������ĳ����������꣬����������ڸ÷����ϵ�����
int Position2Index(float Position) {
	return (Position - MIN) / cube_l;
}


//����������
float Distance(float a[3], float b[3])
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
		(a[1] - b[1]) * (a[1] - b[1]) +
		(a[2] - b[2]) * (a[2] - b[2]));
}

//��ԭ�����
float Distance_0(float x, float y, float z) {
	return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

//�����㣬����ֱ�ߵĲ�������
//	(x-x0)/m = (y-y0)/n = (z-z0)/p
//�˴�ȡm=1
void GenerateLine(float*position1, float*position2, float*paraments) {

	float delta_x = position1[0] - position2[0];
	float delta_y = position1[1] - position2[1];
	float delta_z = position1[2] - position2[2];

	paraments[0] = 1;
	paraments[1] = delta_y / delta_x;
	paraments[2] = delta_z / delta_x;
}

//����һ��������start���������һ��������end
void GeneratePoint(float* start, float* end, float* paraments) {



}

//end�Ƿ񳬳���Χ�з�Χ
bool Inbox(float* end) {

}


//��������ֵ
void GenerateRGBA(float* origin, float* RGBA, float* paraments) {
	float start[3] = { origin[0],origin[1],origin[2] };
	float end[3] = {0,0,0};

	float c_in[3] = { color_sphere[0], color_sphere[1] ,color_sphere[2] };
	float a_in = a_sphere;

	float c_now[3] = {0,0,0};
	float a_now = 0;

	// ���μ����������
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



//������Χ�У���������
void calculate_voxel() {
	float x = MIN;
	float y = MIN;
	float z = MIN;
	//��������
	for (int i = 0; i < NUM; i++, x += cube_l) {	//ѭ��NUM��
		y = MIN;
		for (int j = 0; j < NUM; j++, y += cube_l) {
			z = MIN;
			for (int k = 0; k < NUM; k++, z += cube_l) {
				array[i][j][k].x = x;	//����
				array[i][j][k].y = y;
				array[i][j][k].z = z;
				array[i][j][k].distance = Distance_0(x, y, z);	//�����ľ���
				array[i][j][k].in_out = INOUT(array[i][j][k].distance);	//���ڻ�������
			}
		}
	}
}

//��ǰ����Ⱦ����,������ɫ�Ͳ�͸����
void RayCasting() {

	int count = 0;

	//����ȡֵ
	float RGBA[4] = { 0,0,0,0 };

	//ֱ�߲���
	float paraments[3] = { 1,0,0 };

	//��������
	float image_position[3] = { image[0],SIZE*WIDTH/2,SIZE*HEIGH/2};

	//�������ص�
	for (int i = 0; i < HEIGH; i++) {
		for (int j = 0; j < WIDTH; j++) {
			//����ֱ�߲�������
			GenerateLine(eye, image_position, paraments);

			//��������ֵ
			GenerateRGBA(image_position, RGBA, paraments);

			//��������
			Image[count+0] = RGBA[0];
			Image[count+1] = RGBA[1];
			Image[count+2] = RGBA[2];
			Image[count+3] = RGBA[3];
			count += 4;

			//��һ��������
			image_position[1] += SIZE;
			image_position[2] += SIZE;
		}
	}
}




//ʹ��openglչʾ
void display_voxel() {

	glClear(GL_COLOR_BUFFER_BIT);//����ɫ������
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
	glutCreateWindow("��λ��");
	glutDisplayFunc(display_voxel);
	glutMainLoop();
	int i = 0;
	scanf("%d", &i);
	scanf("%d", &i);

	return 0;
}