#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

const int width = 800;
const int height = 800;
Model* model = NULL;
vec3 light_dir{ 0,0,-1 };//默认光源
TGAImage texture(1024, 1024, TGAImage::RGB);//漫反射贴图

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

void line(const vec2& p1, const vec2& p2, TGAImage& image, TGAColor color)
{
    Bresenham(p1.x, p1.y, p2.x, p2.y, image, color);
}

vec3 barycentric(vec3 A, vec3 B, vec3 C, vec3 P) {
    vec3 s[2];
    for (int i = 2; i--; ) {
        s[i][0] = C[i] - A[i];
        s[i][1] = B[i] - A[i];
        s[i][2] = A[i] - P[i];
    }
    vec3 u = cross(s[0], s[1]);
    if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return vec<3>{1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z};
    return vec<3>{-1, 1, 1}; // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(vec3* pts, vec2* uvs, float* zbuffer, TGAImage& image, float intensity) {
    vec<2> bboxmin{ std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
    vec<2> bboxmax{ -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() };
    vec<2> clamp{ image.get_width() - 1, image.get_height() - 1 };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j])>0.f? std::min(bboxmin[j], pts[i][j]):0.f;
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    vec3 P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
            vec3 bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

            P.z = 0;
            for (int i = 0; i < 3; i++) 
                P.z += pts[i][2] * bc_screen[i];

            if (zbuffer[int(P.x + P.y * width)] < P.z) {
                zbuffer[int(P.x + P.y * width)] = P.z;
                float u = 0, v = 0;
                for (int i = 0; i < 3; i++) {
                    u += bc_screen[i] * uvs[i].x;
                    v += bc_screen[i] * uvs[i].y;
                }
                TGAColor texColor = texture.get(u * texture.get_width(), (1-v) * texture.get_height());
                image.set(P.x, P.y, TGAColor(texColor.r * intensity, texColor.g * intensity, texColor.b * intensity, 255));
            }
        }
    }
}

vec3 world2screen(vec3 v) {
    vec3 out;
    out.x = int((v.x + 1.) * width / 2. + .5);
    out.y = int((v.y + 1.) * height / 2. + .5);
    out.z = v.z;
    return out;
}

//Row scan 
void triangle_row_scan(vec2 t0, vec2 t1, vec2 t2, TGAImage& image, TGAColor color) {
    //冒泡排序,实现 t0.y < t1.y < t2.y
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y; 			//总高度
    for (int y = t0.y; y <= t1.y; y++) { 		//下半部分绘制
        int segment_height = t1.y - t0.y + 1;	//t0到t1的距离包括t1这一行 
        float alpha = (float)(y - t0.y) / total_height; 	// alpha为y轴方向上，绘制的第y行占总高度的比例
        float beta = (float)(y - t0.y) / segment_height; 	// 当心除0错，beta为y轴方向上，绘制的第y行占t1到t的高度的比例
        vec2 A = t0 + (t2 - t0) * alpha;
        vec2 B = t0 + (t1 - t0) * beta;
        if (A.x > B.x) std::swap(A, B); 		 	//从左到右绘制
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// 注意，由于int强制转换t0.y+i！=A.y
        }
    }
    for (int y = t1.y; y <= t2.y; y++) { 		//上半部分绘制
        int segment_height = t2.y - t1.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t1.y) / segment_height; // 当心除0错
        vec2 A = t0 + (t2 - t0) * alpha;
        vec2 B = t1 + (t2 - t1) * beta;
        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++) {
            image.set(j, y, color); 			// 注意，由于int强制转换t0.y+i！=A.y
        }
    }
}

int main(int argc, char** argv) {
    const char* filename = "obj/african_head_diffuse.tga";
    texture.read_tga_file(filename);
    

    TGAImage image(width, height, TGAImage::RGB);
    model = new Model("obj/african_head.obj");
    
    float* zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());//初始化zbuffer每个像素为负无穷

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        std::vector<int> texIndex = model->uv_indices(i);
        vec3 world_coords[3];//三角形的世界坐标
        vec3 pts[3];
        vec2 uvs[3];
        for (int j = 0; j < 3; j++) {
            vec3 v = model->vert(face[j]);
            vec2 uv = model->uv(texIndex[j]);
            world_coords[j] = v;
            pts[j] = world2screen(v);
            uvs[j] = uv;
        }
        vec3 n = cross(world_coords[2] - world_coords[0] , world_coords[1] - world_coords[0]);//三角形法线
        float intensity = n.normalized() * light_dir;//向量单位化,点乘光照得到每个面接受到的光照角度比例
        if (intensity > 0) {//大于0显示，小于0不显示，并根据点乘出来的强度结果影响颜色
            triangle(pts, uvs, zbuffer, image, intensity);
        }
    }

    image.flip_vertically();
    image.write_tga_file("Lec03_zbuffer_diffuse_light.tga");
    return 0;
}