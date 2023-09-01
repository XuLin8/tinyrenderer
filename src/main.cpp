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
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) { // 如果斜率大于1，我们交换 x 和 y
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) { // 如果起点在终点的右边，我们交换起点和终点
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2; // 用于判断是否需要增加 y 的错误值
    int error2 = 0; // 当前的错误值
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color); // 如果斜率大于1，交换 x 和 y 的坐标
        }
        else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) { // 如果错误累积超过 dx，需要增加 y 值
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
    //冒泡排序,实现 t0.y < t1.y < t2.y
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y; 			//总高度
    for (int y = t0.y; y <= t1.y; y++) { 		//下半部分绘制
        int segment_height = t1.y - t0.y + 1;	//t0到t1的距离包括t1这一行 
        float alpha = (float)(y - t0.y) / total_height; 	// alpha为y轴方向上，绘制的第y行占总高度的比例
        float beta = (float)(y - t0.y) / segment_height; 	// 当心除0错，beta为y轴方向上，绘制的第y行占t1到t的高度的比例
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t0 + (t1 - t0) * beta;
        if (A.x > B.x) std::swap(A, B); 		 	//从左到右绘制
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// 注意，由于int强制转换t0.y+i！=A.y
        }
    }
    for (int y = t1.y; y <= t2.y; y++) { 		//上半部分绘制
        int segment_height = t2.y - t1.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t1.y) / segment_height; // 当心除0错
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t1 + (t2 - t1) * beta;
        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// 注意，由于int强制转换t0.y+i！=A.y
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