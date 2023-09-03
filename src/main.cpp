#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

const int width = 800;
const int height = 800;
Model* model = NULL;

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

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    for (int i = 2; i--; ) {
        s[i][0] = C[i] - A[i];
        s[i][1] = B[i] - A[i];
        s[i][2] = A[i] - P[i];
    }
    Vec3f u = s[0] ^ s[1];
    if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(Vec3f* pts, float* zbuffer, TGAImage& image, TGAColor color) {
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
            Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            P.z = 0;
            for (int i = 0; i < 3; i++) P.z += pts[i][2] * bc_screen[i];
            if (zbuffer[int(P.x + P.y * width)] < P.z) {
                zbuffer[int(P.x + P.y * width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
}

//Row scan 
void triangle_row_scan(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
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

    Vec3f light_dir{ 0,0,-1 };//Ĭ�Ϲ�Դ

    TGAImage image(width, height, TGAImage::RGB);
    model = new Model("../obj/african_head.obj");
    float* zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());//��ʼ��zbufferÿ������Ϊ������

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3f world_coords[3];//�����ε���������
        Vec3f pts[3];
        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            world_coords[j] = v;
            pts[j] = world2screen(v);
        }
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);//�����η���
        n.normalize();//������λ��
        float intensity = n * light_dir;//��˹��յõ�ÿ������ܵ��Ĺ��սǶȱ���
        if (intensity > 0) {//����0��ʾ��С��0����ʾ�������ݵ�˳�����ǿ�Ƚ��Ӱ����ɫ
            triangle(pts, zbuffer, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
        }
    }

    image.flip_vertically();
    image.write_tga_file("Lec03_zbuffer.tga");
    return 0;
}