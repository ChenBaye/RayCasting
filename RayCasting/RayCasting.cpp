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
#define DIS		0.01		//��������
#define WIDTH	400	//��Ļ���
#define HEIGH	400	//��Ļ�߶�
#define SIZE		0.0025		//�������صĴ�С
#define INOUT(distance) (Radius-distance)>EPSILON?-1:(distance-Radius)>EPSILON?1:0	
//������distance�ĵ��Ƿ�������	-1-�� 1-�� 0-����

static float center[3] = { 0,0,0 };	//��������

static float eye[3] = { -3,0,0 };	//�ӵ�����

static float image[3] = {-2,0,0};	//��Ļ�������꣬��x�ᴹֱ

static float color_void[3] = { 0,0,0 };	//���������ɫ

static float a_void = 0.005;	//�����ⲻ͸����

static float color_sphere[3] = { 1,0,0 };	//�������ɫ

static float a_sphere = 0.015;	//���岻͸����


//ÿ��ǰ��DIS��xyz������
float delta_x;
float delta_y;
float delta_z;


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


//�����Բ�ֵ
float Linear_Interpolation(float value_0, float position_0, float value_1, float position_1, float position) {

	return value_0 + (value_1 - value_0) * (position - position_0) / (position_1 - position_0);

}

//�����Բ�ֵ����ɫ,���˵�Ϊ���أ�����Ϊx����
void Linear_Interpolation_C(Voxel v1, Voxel v2, float position, float*c, char option='x') {
	if (option == 'x') {
		for (int i = 0; i < 3; i++) {
			c[i] = Linear_Interpolation(
				v1.in_out == -1 ? color_sphere[i] : color_void[i],
				v1.x,
				v2.in_out == -1 ? color_sphere[i] : color_void[i],
				v2.x,
				position
			);
		}
	}
}


//�����Բ�ֵ��͸����,���˵�Ϊ���أ�����Ϊx����
void Linear_Interpolation_A(Voxel v1, Voxel v2, float position, float* a, char option = 'x') {
	if (option == 'x') {
		
		*a = Linear_Interpolation(
			v1.in_out == -1 ? a_sphere : a_void,
			v1.x,
			v2.in_out == -1 ? a_sphere : a_void,
			v2.x,
			position
		);

		/*if (*a < 0) {
			printf("%f %f\n", v1.x, v2.x);
		}*/
		
	}
}


//�Ƿ���Ҫ��ֵ
int NeedInterpolation(Voxel*Voxel_8) {
	int in = 0;
	int out = 0;
	int temp = 0;
	for (int i = 0; i < 8; i++) {
		temp = Voxel_8[i].in_out;
		if (temp == -1) {

			in = 1;
		}
		else if (temp == 1) {
			out = 1;
		}
	}
	
	if (in != 0 && out != 0) {
		return 0;	//��Ҫ��ֵ
	}
	else if (in != 0) {
		//printf("in\n");
		return -1;	//ȫ��λ�����ڲ�(������)
	}	
	else {			
		return 1;	//ȫ��λ�����ⲿ(������)
	}
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

	float sub_x = position2[0] - position1[0];
	float sub_y = position2[1] - position1[1];
	float sub_z = position2[2] - position1[2];

	paraments[0] = 1;
	paraments[1] = sub_y / sub_x;
	paraments[2] = sub_z / sub_x;

	//����x,y,z����
	float SUM_2 = paraments[0] * paraments[0] + paraments[1] * paraments[1] + paraments[2] * paraments[2];
	float temp = sqrt(SUM_2);

	delta_x = DIS * paraments[0] / temp;
	delta_y = DIS * paraments[1] / temp;
	delta_z = DIS * paraments[2] / temp;
	//printf("%f %f %f\n", delta_x, delta_y, delta_z);
	//printf("GenerateLine over\n");
}

//����һ��������start���������һ��������end
void GeneratePoint(float* start, float* end, float* paraments) {
	//	(x-x0)/m = (y-y0)/n = (z-z0)/p
	//	�˴�ȡm=1
	
	end[0] = start[0] + delta_x;
	end[1] = start[1] + delta_y;
	end[2] = start[2] + delta_z;
}

//end�Ƿ񳬳���Χ�з�Χ
bool Inbox(float* end) {
	for (int i = 0; i < 3; i++) {
		if (end[i]<MIN || end[i]>MAX) {
			return false;
		}
	}

	
	return true;
}

//λ�ڰ�Χ����image֮��
bool Beforebox(float* end) {
	if (end[0] >= image[0] && end[0] < MIN) {
		return true;
	}
	else {
		return false;
	}
}

//������һ��������
void GenerateStart(float *start) {
	
	while (true) {
		if(Inbox(start)){
			//printf("inbox\n");
			break;
		}
		else if (Beforebox(start)) {
			start[0] += delta_x;
			start[1] += delta_y;
			start[2] += delta_z;
			//printf("beforebox\n");
		}
		else {	//������Χ��
			//printf("outbox\n");
			break;
		}
	}
}

// ����������Ԫ����ĳ����������С����
float GetMin(float x) {
	return ((int)(x / cube_l)-1)*cube_l;
}

// ����������Ԫ����ĳ���������������
float GetMax(float x) {
	return ((int)(x / cube_l))*cube_l;
}

//����c_now��a_now
void Calculate_C_A(float*end, float*c_now, float*a_now) {
	//�����8����������
	float x_min = GetMin(end[0]);
	float x_max = GetMax(end[0]);

	float y_min = GetMin(end[1]);
	float y_max = GetMax(end[1]);

	float z_min = GetMin(end[2]);
	float z_max = GetMax(end[2]);

	//�洢end��������Ԫ�İ˸�����
	Voxel Voxel_8[8] = {
		array[Position2Index(x_min)][Position2Index(y_min)][Position2Index(z_min)],
		array[Position2Index(x_max)][Position2Index(y_min)][Position2Index(z_min)],
		array[Position2Index(x_max)][Position2Index(y_max)][Position2Index(z_min)],
		array[Position2Index(x_min)][Position2Index(y_max)][Position2Index(z_min)],
		array[Position2Index(x_min)][Position2Index(y_min)][Position2Index(z_max)],
		array[Position2Index(x_max)][Position2Index(y_min)][Position2Index(z_max)],
		array[Position2Index(x_max)][Position2Index(y_max)][Position2Index(z_max)],
		array[Position2Index(x_min)][Position2Index(y_max)][Position2Index(z_max)],
	};
	//printf("%d %d %d\n", Position2Index(x_min), Position2Index(y_min), Position2Index(z_min));
	
	//�Ƿ���Ҫ��ֵ
	int result = NeedInterpolation(Voxel_8);
	switch (result) {
		case -1:		//ȫ���ڲ�
			c_now[0] = color_sphere[0];
			c_now[1] = color_sphere[1];
			c_now[2] = color_sphere[2];
			*a_now = a_sphere;
			break;
		case 1:		//ȫ���ⲿ
			c_now[0] = color_void[0];
			c_now[1] = color_void[1];
			c_now[2] = color_void[2];
			*a_now = a_void;
			break;
		case 0:		//�������⣬��Ҫ��ֵ
			//�������Բ�ֵ��Ϊ�ߴε����Բ�ֵ
			//��ȡ��end�㴹ֱ��x���ƽ���ȡ��Ԫ���ĸ���
			float c_p1[3], a_p1;
			float c_p2[3], a_p2;
			float c_p3[3], a_p3;
			float c_p4[3], a_p4;


			Linear_Interpolation_C(Voxel_8[0], Voxel_8[1], end[0], c_p1, 'x');
			Linear_Interpolation_C(Voxel_8[2], Voxel_8[3], end[0], c_p2, 'x');
			Linear_Interpolation_C(Voxel_8[6], Voxel_8[7], end[0], c_p3, 'x');
			Linear_Interpolation_C(Voxel_8[4], Voxel_8[5], end[0], c_p4, 'x');
			
			
			Linear_Interpolation_A(Voxel_8[0], Voxel_8[1], end[0], &a_p1, 'x');
			Linear_Interpolation_A(Voxel_8[2], Voxel_8[3], end[0], &a_p2, 'x');
			Linear_Interpolation_A(Voxel_8[6], Voxel_8[7], end[0], &a_p3, 'x');
			Linear_Interpolation_A(Voxel_8[4], Voxel_8[5], end[0], &a_p4, 'x');

			//�ٽ������ε����Բ�ֵ
			float c_p5[3], a_p5;
			float c_p6[3], a_p6;

			for (int i = 0; i < 3; i++) {
				c_p5[i] = Linear_Interpolation(c_p1[i], z_min, c_p4[i], z_max, end[2]);
				c_p6[i] = Linear_Interpolation(c_p2[i], z_min, c_p3[i], z_max, end[2]);
			}
			a_p5 = Linear_Interpolation(a_p1, z_min, a_p4, z_max, end[2]);
			a_p6 = Linear_Interpolation(a_p2, z_min, a_p3, z_max, end[2]);

			//���һ�ε����Բ�ֵ
			for (int i = 0; i < 3; i++) {
				c_now[i] = Linear_Interpolation(c_p5[i], y_min, c_p6[i], y_max, end[1]);
			}
			*a_now = Linear_Interpolation(a_p5, y_min, a_p6, y_max, end[1]);

			//printf("%f\n", *a_now);
			break;
	}

}


//��������ֵ
void GenerateRGBA(float* origin, float* RGBA, float* paraments) {
	//printf("into GenerateRGBA\n");
	//printf("%f %f %f\n", origin[0], origin[1], origin[2]);
	float start[3] = { origin[0],origin[1],origin[2] };
	float end[3] = {0,0,0};

	//����������ǰ
	float c_in[3] = { 0, 0 ,0 };
	float a_in = 0;

	//������
	float c_now[3] = {0,0,0};
	float a_now = 0;

	//�����������
	float c_out[3] = { 0,0,0 };
	float a_out = 0;

	//test
	start[0] = -2;
	start[1] = 0;
	start[2] = 0;
	delta_x = DIS;
	delta_y = 0;
	delta_z = 0;
	//test


	//�����һ��������
	GenerateStart(start);
	//printf("%f %f %f\n", start[0], start[1], start[2]);

	while (true) {
		// �����������

		GeneratePoint(start, end, paraments);
		//printf("%f %f %f\n", end[0], end[1], end[2]);

		if (Inbox(end) && a_out<1) {		//���ܳ�����Χ�У����ۼƲ�͸���Ȳ��ܳ���1
			Calculate_C_A(end, c_now, &a_now);	//����c_now��a_now
			printf("%f %f %f\n", end[0], end[1], end[2]);
			printf("%f %f %f %f\n", c_now[0], c_now[1], c_now[2], a_now);

			a_out = a_in + a_now*(1-a_in);	//��͸����A
			c_out[0] = c_in[0]*a_in + c_now[0]*a_now*(1 - a_in);	//��ɫR
			c_out[1] = c_in[1]*a_in + c_now[1]*a_now*(1 - a_in);	//��ɫG
			c_out[2] = c_in[2]*a_in + c_now[2]*a_now*(1 - a_in);	//��ɫB

			//printf("%f %f %f %f\n", c_out[0], c_out[1], c_out[2], a_out);

			a_in = a_out;
			c_in[0] = c_out[0];
			c_in[1] = c_out[1];
			c_in[2] = c_out[2];

			start[0] = end[0];
			start[1] = end[1];
			start[2] = end[2];
		}
		else {
			break;
		}
	}

	RGBA[0] = c_out[0];
	RGBA[1] = c_out[1];
	RGBA[2] = c_out[2];
	RGBA[3] = a_out;
	//printf("%f %f %f %f\n", c_out[0], c_out[1], c_out[2], a_out);

	//printf("GenerateRGBA over\n");
}



//������Χ�У���������
void calculate_voxel() {
	float x = MIN;
	float y = MIN;
	float z = MIN;
	//��������
	for (int i = 0; i < NUM; i++, z += cube_l) {	//ѭ��NUM��
		y = MIN;
		for (int j = 0; j < NUM; j++, y += cube_l) {
			x = MIN;
			for (int k = 0; k < NUM; k++, x += cube_l) {
				array[i][j][k].x = x;	//����
				array[i][j][k].y = y;
				array[i][j][k].z = z;
				array[i][j][k].distance = Distance_0(x, y, z);	//�����ľ���
				array[i][j][k].in_out = INOUT(array[i][j][k].distance);	//���ڻ�������
			}
		}
	}

	//printf("calculate_voxel over\n");
}

//��ǰ����Ⱦ����,������ɫ�Ͳ�͸����
void RayCasting() {

	int count = 0;

	//����ȡֵ
	float RGBA[4] = { 0,0,0,0 };

	//ֱ�߲���
	float paraments[3] = { 1,0,0 };

	//��������
	float image_position[3] = { image[0],-SIZE*WIDTH/2,-SIZE*HEIGH/2};

	//�������ص�
	for (int i = 0; i < HEIGH; i++) {
		printf("%d\n",i);
		image_position[1] = -SIZE * WIDTH / 2;
		for (int j = 0; j < WIDTH; j++) {
			//printf("%d\n", j);
			//����ֱ�߲�������
			GenerateLine(eye, image_position, paraments);
			//printf("%f %f %f\n", image_position[0], image_position[1], image_position[2]);

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
		}
		image_position[2] += SIZE;
	}

	printf("RayCasting over\n");
}




//ʹ��openglչʾ
void display_voxel() {


	
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(WIDTH, HEIGH, GL_RGBA, GL_FLOAT, Image);//ʹ��OpenGL�Ļ�ͼ����
	glFlush();
	
}


int main(int argc, char** argv) {

	calculate_voxel();
	RayCasting();



	printf("%f", Image[200]);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
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