#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const int width = 200;
const int height = 200;

//Bresenham
void Bresenham(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) { // ���б�ʴ���1�����ǽ��� x �� y
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) { // ���������յ���ұߣ����ǽ��������յ�
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2; // �����ж��Ƿ���Ҫ���� y �Ĵ���ֵ
    int error2 = 0; // ��ǰ�Ĵ���ֵ
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color); // ���б�ʴ���1������ x �� y ������
        }
        else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) { // ��������ۻ����� dx����Ҫ���� y ֵ
            y += (y1 > y0 ? 1 : -1);
            error2 -= dx * 2;
        }
    }
}

void line(const Vec2i& p1, const Vec2i& p2, TGAImage& image, TGAColor color)
{
    Bresenham(p1.x, p1.y, p2.x, p2.y, image, color);
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    //ð������,ʵ�� t0.y < t1.y < t2.y
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y; 			//�ܸ߶�
    for (int y = t0.y; y <= t1.y; y++) { 		//�°벿�ֻ���
        int segment_height = t1.y - t0.y + 1;	//t0��t1�ľ������t1��һ�� 
        float alpha = (float)(y - t0.y) / total_height; 	// alphaΪy�᷽���ϣ����Ƶĵ�y��ռ�ܸ߶ȵı���
        float beta = (float)(y - t0.y) / segment_height; 	// ���ĳ�0��betaΪy�᷽���ϣ����Ƶĵ�y��ռt1��t�ĸ߶ȵı���
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t0 + (t1 - t0) * beta;
        if (A.x > B.x) std::swap(A, B); 		 	//�����һ���
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// ע�⣬����intǿ��ת��t0.y+i��=A.y
        }
    }
    for (int y = t1.y; y <= t2.y; y++) { 		//�ϰ벿�ֻ���
        int segment_height = t2.y - t1.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t1.y) / segment_height; // ���ĳ�0��
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t1 + (t2 - t1) * beta;
        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// ע�⣬����intǿ��ת��t0.y+i��=A.y
        }
    }
}

int main(int argc, char** argv) {
    
    TGAImage image(width, height, TGAImage::RGB);

    Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
    Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
    Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };

    triangle(t0[0], t0[1], t0[2], image, red);
    triangle(t1[0], t1[1], t1[2], image, white);
    triangle(t2[0], t2[1], t2[2], image, green);

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("triangle.tga");
    return 0;
}