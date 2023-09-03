# 指定 CMake 的最低版本要求
cmake_minimum_required(VERSION 3.0)

# 设置 C++ 标准为 C++17
set(CMAKE_CXX_STANDARD 17)

# 禁用 ALL_BUILD 和 ZERO_CHECK 目标
set(CMAKE_SUPPRESS_REGENERATION true)

# 为项目命名
project(tinyrenderer)

# 将 tinyrenderer 设置为默认启动项目
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT tinyrenderer)

# 查找并添加所有的源文件
file(GLOB SRCS "src/*.cpp" "src/*.h")

# 创建可执行文件
add_executable(tinyrenderer ${SRCS})

# 添加头文件搜索路径
target_include_directories(tinyrenderer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # 添加头文件目录
)

# 输出一些构建信息
message("C++ 编译器: ${CMAKE_CXX_COMPILER_ID}")
message("C++ 标准: C++${CMAKE_CXX_STANDARD}")

