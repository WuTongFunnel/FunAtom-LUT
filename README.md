## FunAtom-LUT

A multi-threaded C++ 3D LUT baker for physically based atmospheric scattering.

This program performs raymarching integration with recursive multiple scattering, earth horizon occlusion, and Hammersley spherical sample for indirect sky light.
It precomputes a 128³ RGBAHalf volume lookup table (`atom.bin`) storing sky radiance and transmittance data.
The baked binary LUT can be imported into real-time renderers, including Minecraft Iris shaders, GLSL/HLSL shader pipelines and custom soft renderers, to accelerate sky & atmosphere rendering at runtime.

### Features

- Raymarching atmospheric scattering solver with recursive multi-scatter
- Earth occlusion & sun disk solid-angle visibility calculation
- Hammersley low-discrepancy spherical sampling for indirect light
- Multi-threaded baking with atomic progress counter
- Built-in float32 to float16 (half-precision) conversion
- Output binary: 128³ RGBAHalf `atom.bin`

### Build

Requires C++17 compliant compiler.
Dependencies: Standard library only, no external third-party libraries.


> 
> Note: `#include <graphics.h>` is only retained for your local soft renderer prototype; it is unused in the LUT baking logic and can be safely removed when compiling on Codespaces / Linux.

# FunAtom-LUT

基于物理大气散射的 C++ 多线程 3D LUT 烘焙器。

程序使用光线步进积分实现递归多次散射，支持地球地平线遮挡、Hammersley 球面低差异采样求解天空间接光。
预计算生成 128³ RGBAHalf 体积查找表 `atom.bin`，存储天空辐亮度与透射率数据。
烘焙完成的二进制 LUT 可接入实时渲染管线，例如 Minecraft Iris 着色器、GLSL/HLSL 着色器以及自研软渲染器，用于在运行时加速大气天空渲染。

### 特性

- 支持递归多次散射的光线步进大气散射求解器
- 地球遮挡与太阳圆盘立体角可见性计算
- Hammersley 低差异球面采样用于间接光照
- 多线程烘焙，原子计数器实时输出进度
- 内置 float32 转 float16 半精度转换
- 输出二进制文件：128³ RGBAHalf `atom.bin`

### 编译

需要支持 C++17 的编译器。
仅依赖 C++ 标准库，无第三方依赖。
