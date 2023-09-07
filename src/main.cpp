#include <vector>
#include <iostream>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model* model = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir(1, 1, 1);
Vec3f       eye(0, -1, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

// 自定义着色器，继承自IShader
struct GouraudShader : public IShader {
    Vec3f varying_intensity; // 由顶点着色器写入，片段着色器读取

    // 顶点着色器
    virtual Vec4f vertex(int iface, int nthvert) {
        // 从.obj文件中读取顶点坐标
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
        // 将顶点坐标变换到屏幕坐标
        gl_Vertex = Viewport * Projection * ModelView * gl_Vertex;
        // 计算漫反射光照强度
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * light_dir);
        return gl_Vertex;
    }

    // 片段着色器
    virtual bool fragment(Vec3f bar, TGAColor& color) {
        // 插值计算当前像素的光照强度
        float intensity = varying_intensity * bar;
        // 根据光照强度设置颜色
        color = TGAColor(255, 255, 255) * intensity;
        // 不丢弃当前像素
        return false;
    }
};

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head.obj");
    }

    // 设置相机参数
    lookat(eye, center, up);
    // 设置视口参数
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    // 设置透视投影
    projection(-1.f / (eye - center).norm());
    light_dir.normalize();

    // 创建图像和深度缓冲
    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    // 创建着色器
    GouraudShader shader;
    for (int i = 0; i < model->nfaces(); i++) {
        Vec4f screen_coords[3];
        for (int j = 0; j < 3; j++) {
            screen_coords[j] = shader.vertex(i, j);
        }
        // 绘制三角形
        triangle(screen_coords, shader, image, zbuffer);
    }

    // 图像翻转以将原点放在图像的左下角
    image.flip_vertically();
    zbuffer.flip_vertically();
    // 保存图像和深度缓冲
    image.write_tga_file("imgs/Lec06_output.tga");
    zbuffer.write_tga_file("imgs/Lec06_zbuffer.tga");

    delete model;
    return 0;
}
