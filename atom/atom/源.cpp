#include <graphics.h>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <iostream>
using namespace std;

#define PI 3.141592653589793f

struct vec3
{
    float x, y, z;
    vec3() : x(0), y(0), z(0) {}
    vec3(float v) : x(v), y(v), z(v) {}
    vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    vec3 operator+(const vec3& o) const { return vec3(x + o.x, y + o.y, z + o.z); }
    vec3 operator-(const vec3& o) const { return vec3(x - o.x, y - o.y, z - o.z); }
    vec3 operator*(const vec3& o) const { return vec3(x * o.x, y * o.y, z * o.z); }
    vec3 operator*(float s) const { return vec3(x * s, y * s, z * s); }
    vec3 operator/(float s) const { return vec3(x / s, y / s, z / s); }
    // 支持 float * vec3
    friend vec3 operator*(float s, const vec3& v)
    {
        return v * s;
    }
};

inline float dot(const vec3& a, const vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float length(const vec3& v)
{
    return sqrtf(dot(v, v));
}
inline vec3 normalize(const vec3& v)
{
    float l = length(v);
    if (l < 1e-6f) return vec3(0, 0, 0);
    return v / l;
}
inline float clamp(float v, float lo, float hi)
{
    return max(lo, min(hi, v));
}
inline float radians(float deg)
{
    return deg * PI / 180.0f;
}

// ===================== 全局参数（对应GLSL uniform） =====================
const float mie_g = 0.375f;
const float atom_height = 100000.0f;
const float earth_factor = 1.0f;
const float earth_r = 6371000.0f;
const float MC_r = earth_r / earth_factor;
const float MC_atom_height = atom_height / earth_factor;
const float MC_ar = MC_r + MC_atom_height;
const float MC_sealevel = 63.0f;

float sun_atten_r = 0.955f;
float sun_atten_g = 0.897f;
float sun_atten_b = 0.768f;
vec3 cameraPosition=vec3(0,MC_sealevel,0);
vec3 E=vec3(0, MC_sealevel - cameraPosition.y - MC_r, 0);
float sunAzimuth = 0.0f;
float sunElevation = 0.4f;
const float sunStep = 0.02f;

// ===================== 函数声明 =====================
vec3 scatterf(float l);
float intersectEarthT1(vec3 p, vec3 beginpostion);
float intersectEarthT2(vec3 p,vec3 beginpostion);
float sunRayAtmLength(vec3 P, vec3 sunDir);
float capSolidAngle(float halfAngle);
float intersectCapsSolidAngle(float alpha, float beta, float gamma);
float sunCanReachPoint(vec3 P, vec3 sunDir);
vec3 HammersleySphereDir(uint32_t i, uint32_t sampleTotal);

// ===================== 函数实现 =====================
float fog_k(vec3 sun_world)
{
    float sun_theta_c = dot(sun_world, vec3(0, 1, 0));
    float sun_theta_s = sqrt(1 - pow(sun_theta_c, 2));
    float fog_world_fun_max = 1.05f;
    float fog_world_fun_min = 1.0f;
    float fog_world_fun_speed = 4.0f;
    float fog_world_fun = (fog_world_fun_max - fog_world_fun_min) * powf(sun_theta_s, 2.0f * fog_world_fun_speed) + fog_world_fun_min;
    return fog_world_fun;
}

float fog_kk(float l2,vec3 sun_world)
{
    float mix_sky_k = 1.0f - powf(fog_k(sun_world), -l2);
    return mix_sky_k;
}

vec3 scatterf(float l)
{
    return vec3(powf(sun_atten_r, l), powf(sun_atten_g, l), powf(sun_atten_b, l));
}

float intersectEarthT1(vec3 p,vec3 beginpostion)
{
 
    vec3 oc = beginpostion - E;
    float b = 2.0f * dot(p, oc);
    float c = dot(oc, oc) - MC_ar * MC_ar;
    float delta = b * b - 4.0f * c;
    if (delta < 0.0f) return -1.0f;
    float sqrtDelta = sqrtf(delta);
    float t1 = (-b - sqrtDelta) / 2.0f;
    return max(0.0f, t1);
}

float intersectEarthT2(vec3 p, vec3 beginpostion)
{
    vec3 oc = beginpostion - E;
    float b = 2.0f * dot(p, oc);
    float c = dot(oc, oc) - MC_ar * MC_ar;
    float delta = b * b - 4.0f * c;
    if (delta < 0.0f) return -1.0f;
    float sqrtDelta = sqrtf(delta);
    float t2 = (-b + sqrtDelta) / 2.0f;
    return max(0.0f, t2);
}

float sunRayAtmLength(vec3 P, vec3 sunDir)
{
    vec3 rel = P - E;
    float A = dot(sunDir, rel);
    float lenRelSq = dot(rel, rel);
    float R = MC_ar;
    float delta = A * A - (lenRelSq - R * R);
    if (delta < 0.0f) return 0.0f;
    float sqrtDelta = sqrtf(delta);
    float t = -A - sqrtDelta;
    return max(max(-A - sqrtDelta, -A + sqrtDelta), 0.0f);
}

const float EPS = 1e-6f;
float capSolidAngle(float halfAngle)
{
    return 2.0f * PI * (1.0f - cosf(halfAngle));
}

float intersectCapsSolidAngle(float alpha, float beta, float gamma)
{
    alpha = clamp(alpha, 0.0f, PI);
    beta = clamp(beta, 0.0f, PI);
    gamma = clamp(gamma, 0.0f, PI);
    if (alpha >= PI - EPS) return capSolidAngle(beta);
    if (beta >= PI - EPS) return capSolidAngle(alpha);
    if (gamma >= PI - EPS) return 0.0f;
    if (gamma >= alpha + beta - EPS) return 0.0f;
    if (gamma <= fabsf(alpha - beta) + EPS)
        return capSolidAngle(min(alpha, beta));

    float sa = sinf(alpha);
    float sb = sinf(beta);
    float sg = sinf(gamma);
    if (sa < EPS || sb < EPS || sg < EPS) return 0.0f;

    float lambdaA = acosf(clamp(
        (cosf(beta) - cosf(alpha) * cosf(gamma)) / (sa * sg),
        -1.0f, 1.0f));
    float lambdaB = acosf(clamp(
        (cosf(alpha) - cosf(beta) * cosf(gamma)) / (sb * sg),
        -1.0f, 1.0f));
    float lambdaP = acosf(clamp(
        (cosf(gamma) - cosf(alpha) * cosf(beta)) / (sa * sb),
        -1.0f, 1.0f));

    float omega = 2.0f * PI
        - 2.0f * (lambdaP + lambdaA * cosf(alpha) + lambdaB * cosf(beta));
    return clamp(omega, 0.0f, capSolidAngle(min(alpha, beta)));
}

float sunCanReachPoint(vec3 P, vec3 sunDir)
{
    vec3 v = normalize(E - P);
    float earth_distance = length(P - E);
    float earth_halftheta = asinf(MC_r / earth_distance);
    if (earth_distance < MC_r) return 0.0f;
    float sun_halftheta = radians(0.25f);
    float sun_s = capSolidAngle(sun_halftheta);
    float dotVSun = dot(v, sunDir);
    dotVSun = clamp(dotVSun, -1.0f, 1.0f);
    float gamma = acosf(dotVSun);
    float cross_s = intersectCapsSolidAngle(sun_halftheta, earth_halftheta, gamma);
    float sun_visible = sun_s - cross_s;
    return clamp(sun_visible / sun_s, 0.0f, 1.0f);
}

vec3 HammersleySphereDir(uint32_t i, uint32_t sampleTotal)
{
    uint32_t bits = i;
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    float rInv = (float)bits * 2.3283064365386963e-10f;
    float u = (float)i / (float)sampleTotal;
    float v = rInv;
    float phi = 2.0f * PI * u;
    float z = 1.0f - 2.0f * v;
    float r = sqrtf(max(0.0f, 1.0f - z * z));
    return vec3(r * cosf(phi), r * sinf(phi), z);
}

// vec2简易实现
struct vec2 { float x, y; vec2() {} vec2(float a, float b) :x(a), y(b) {} };

int n= 1;
vec3 fogmodel(vec3 pixel_world,vec3 beginpostion,vec3 sun_world,int k)
{
    if (k == n + 1)return vec3(0);
    vec3 allcolor = vec3(0);
   // if (dot(pixel_world, sun_world) > 0.9999)return vec3(1);
    vec3 rm_v = normalize(pixel_world);
    float t1 = intersectEarthT1(rm_v, beginpostion);
    float t2 = intersectEarthT2(rm_v,beginpostion);
float rm_all_l = t2 - t1;
    if (rm_all_l < 0.001f) return allcolor;
    if (rm_all_l == -1.0f) return allcolor;
    float rm_t =128.0f;
    float rm_dx = rm_all_l / rm_t;
    float pixel_sun_theta_c = dot(rm_v, sun_world);
    float rayleigh_phase = (3.0f / (16.0f * PI)) * (1.0f + pixel_sun_theta_c * pixel_sun_theta_c);
    float mie_phase = (1.0f / (4.0f * PI)) * (1.0f - mie_g * mie_g)
        / powf(1.0f + mie_g * mie_g - 2.0f * mie_g * pixel_sun_theta_c, 1.5f);
    float opixel_sun_theta_c = 1.0f;
    float orayleigh_phase = (3.0f / (16.0f * PI)) * (1.0f + opixel_sun_theta_c * opixel_sun_theta_c);
    float omie_phase = (1.0f / (4.0f * PI)) * (1.0f - mie_g * mie_g)
        / powf(1.0f + mie_g * mie_g - 2.0f * mie_g * opixel_sun_theta_c, 1.5f);
    vec3 rm_p = rm_all_l * rm_v +t1 * rm_v + beginpostion;
    vec3 p_skycolor(0, 0, 0);
    for (int i = 1; i <= (int)rm_t; i++)
    {
        if (length(rm_p - E) > MC_ar)
        {
            rm_p = rm_p - rm_v * rm_dx;
            continue;
        }
        vec3 o_color = allcolor * scatterf(rm_dx / MC_atom_height);
        vec3 new_r_color = (allcolor - o_color) * orayleigh_phase;
        o_color = o_color * (1.0f - fog_kk(rm_dx / MC_atom_height,sun_world));
        vec3 new_m_color = o_color * fog_kk(rm_dx / MC_atom_height,sun_world) * omie_phase;
        allcolor = o_color + new_r_color + new_m_color;
        float sun_to_rm_p = sunRayAtmLength(rm_p, sun_world);
            vec3 p_suncolor = vec3(1, 1, 1)* scatterf(sun_to_rm_p / MC_atom_height)
                * (1.0f - fog_kk(sun_to_rm_p / MC_atom_height,sun_world));
                p_suncolor = p_suncolor * sunCanReachPoint(rm_p, sun_world);
            vec3 p_allcolor = p_suncolor;
            vec3 drscateer(
                -logf(sun_atten_r) * rm_dx / MC_atom_height,
                -logf(sun_atten_g) * rm_dx / MC_atom_height,
                -logf(sun_atten_b) * rm_dx / MC_atom_height
            );
            vec3 drcolor = p_allcolor * drscateer * rayleigh_phase;
            float dmscatter = (logf(fog_k(sun_world)) * rm_dx / MC_atom_height);
            vec3 dmcolor = mie_phase * p_suncolor * dmscatter;
            vec3 dcolor = dmcolor + drcolor;
            allcolor = allcolor + dcolor;

            {
                p_skycolor = vec3(0);
                int t = 9;
                if (k < n + 1) {
                    for (int j= 0; j < t; j++)
                    {
                        vec3 p = HammersleySphereDir(j, t);
                        vec3 tempcolor = fogmodel(p, rm_p, sun_world,k + 1);
                        float p_main_theta_c = dot(rm_v, p);
                        float prayleigh_phase = (3 / (16 * PI)) * (1 + p_main_theta_c * p_main_theta_c);
                        float pmie_phase = (1.0 / (4.0 * PI)) * (1.0 - mie_g * mie_g) / pow(1.0 + mie_g * mie_g - 2.0 * mie_g * p_main_theta_c, 1.5);
                        p_skycolor = p_skycolor + tempcolor * drscateer * prayleigh_phase;
                        p_skycolor = p_skycolor + tempcolor * dmscatter * pmie_phase;

                    }
                    p_skycolor = p_skycolor / t;
                    p_skycolor = p_skycolor * 4 * PI;
                }
            }
        allcolor = allcolor+ p_skycolor;
        rm_p = rm_p - rm_v * rm_dx;
    }
    return allcolor;
}

// 存储每个LUT格子结果：rgb=天空辐亮度，a=透射率，你可以按需修改
struct RGBAHalf
{
    uint16_t r, g, b, a;
};

#include <fstream>
#include <cstdint>
#include <vector>
// float32 -> float16 转换
#include <cstdint>

uint16_t float_to_half(float f)
{
    // ========== 在函数内部做数值清洗 ==========
    const float eps = 1e-4f;
    const float halfMax = 65500.0f;

    // 负数归零
    if (f < 0.0f)
        f = 0.0f;
    // 极小值归零
    if (f < eps)
        f = 0.0f;
    // 上限截断，防止Inf
    if (f > halfMax)
        f = halfMax;

    // ========== 下面是你原来的位运算代码，完全不变 ==========
    uint32_t bits = *(uint32_t*)&f;
    uint16_t sign = (bits >> 31) & 0x1;
    uint16_t exp = (bits >> 23) & 0xFF;
    uint32_t mant = bits & 0x007FFFFF;
    uint16_t h_exp, h_mant;
    if (exp == 0)
    {
        // float本身是非正规数，直接输出0，禁止生成half denormal
        h_exp = 0;
        h_mant = 0;
    }
    else if (exp == 255)
    {
        // Inf / NaN
        h_exp = 31;
        h_mant = mant ? 0x200 : 0;
    }
    else
    {
        int e = exp - 127;
        e += 15;
        if (e >= 31)
        {
            // 溢出 → Inf
            h_exp = 31;
            h_mant = 0;
        }
        else if (e <= 0)
        {
            // 指数太小，不再计算denormal，直接归零！重点修改这里
            h_exp = 0;
            h_mant = 0;
        }
        else
        {
            h_exp = static_cast<uint16_t>(e);
            h_mant = static_cast<uint16_t>(mant >> 13);
        }
    }
    return (sign << 15) | (h_exp << 10) | h_mant;
}


#include <thread>
#include <vector>


#include <thread>
#include <vector>
#include <thread>
#include <vector>

// 常量提升到全局
const int LUT_W = 128;
const int LUT_H = 128;
const int LUT_D = 128;
// ========== 在这里调节：每个任务一次性算多少个Z ==========
const int Z_PER_TASK = 4;

#include <atomic>
#include <thread>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>

// 原子计数器，统计已经完成的像素数量
std::atomic<uint64_t> donePixels{ 0 };
// 全局任务：下一组z起始
std::atomic<int> nextZTask{ 0 };

// 单个分片计算任务（你原来的函数，完全不动）
void ComputeZSlice(int zStart, int zEnd, std::vector<RGBAHalf>& outBuffer)
{
    outBuffer.clear();
    int totalPixelInSlice = (zEnd - zStart) * LUT_H * LUT_W;
    outBuffer.reserve(totalPixelInSlice);
    uint64_t localCnt = 0;
    const uint64_t batch =64;
    for (int z = zStart; z < zEnd; z++)
    {
        for (int y = 0; y < LUT_H; y++)
        {
            for (int x = 0; x < LUT_W; x++)
            {
                float u = (float)x / (LUT_W - 1);
                float v = (float)y / (LUT_H - 1);
                float w = (float)z / (LUT_D - 1);
                float mu_v = 2.0f * u - 1.0f;
                float mu_s = 2.0f * v - 1.0f;
                float cosGamma = 2.0f * w - 1.0f;
                float sinThetasSq = max(0.0f, 1.0f - mu_s * mu_s);
                float s_x = sqrtf(sinThetasSq);
                vec3 sun_world = vec3(s_x, mu_s, 0);
                float sinThetavSq = max(0.0f, 1.0f - mu_v * mu_v);
                float sinThetav = sqrtf(sinThetavSq);
                float cosPhi = cosGamma;
                float sinPhi = sqrtf(max(0.0f, 1.0f - cosGamma * cosGamma));
                vec3 view_dir = vec3(
                    sinThetav * cosPhi,
                    mu_v,
                    sinThetav * sinPhi
                );
                vec3 skyColor;
                skyColor = fogmodel(view_dir, vec3(0), sun_world, 1);
                float transmittance = 1.0f;
                RGBAHalf pix;
                pix.r = float_to_half(skyColor.x);
                pix.g = float_to_half(skyColor.y);
                pix.b = float_to_half(skyColor.z);
                pix.a = float_to_half(transmittance);
                outBuffer.push_back(pix);
                localCnt++;
                if (localCnt >= batch)
                {
                    donePixels.fetch_add(batch);
                    localCnt -= batch;
                }
            }
        }
    }
    if (localCnt > 0)
    {
        donePixels.fetch_add(localCnt);
    }
}

int Write3DLUTBin()
{
    // 自动获取可用并发线程数
    unsigned int threadCount = std::thread::hardware_concurrency();
    if (threadCount == 0)
    {
        threadCount = 1;
    }
    std::cout << "自动检测线程数量：" << threadCount << std::endl;
    std::vector<std::thread> workers;

    int zTotal = LUT_D;
    // 存储所有层结果，按z编号保存，保证写入顺序正确
    std::vector<std::vector<RGBAHalf>> lutStorage(zTotal);

    auto workerFunc = [&]()
        {
            while (true)
            {
                // 原子拿任务起始z
                int zStart = nextZTask.fetch_add(Z_PER_TASK);
                if (zStart >= zTotal) break;

                int zEnd = zStart + Z_PER_TASK;
                if (zEnd > zTotal) zEnd = zTotal; // 最后一组不足Z_PER_TASK

                std::vector<RGBAHalf> tmpBuf;
                ComputeZSlice(zStart, zEnd, tmpBuf);

                // 把计算结果拷贝到对应z范围的存储空间
                size_t ptr = 0;
                for (int z = zStart; z < zEnd; z++)
                {
                    lutStorage[z] = std::vector<RGBAHalf>(tmpBuf.begin() + ptr,
                        tmpBuf.begin() + ptr + LUT_W * LUT_H);
                    ptr += LUT_W * LUT_H;
                }
            }
        };

    // 启动所有worker
    for (unsigned int i = 0; i < threadCount; i++)
    {
        workers.emplace_back(workerFunc);
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    uint64_t total = (uint64_t)LUT_W * LUT_H * LUT_D;
    while (true)
    {
        uint64_t finished = donePixels.load();
        if (finished >= total)
            break;
        auto now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(now - startTime).count();
        double progress = static_cast<double>(finished) / total;
        double remain = elapsed / progress - elapsed;
        std::cout << "\rProgress: " << std::fixed << progress * 100
            << "% | Elapsed: " << elapsed << "s | Remain: " << remain << "s";
        std::cout.flush();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    // 等待所有线程执行完毕
    for (auto& t : workers)
    {
        t.join();
    }
    // 主线程按z从小到大写入文件，顺序永远正确
    std::ofstream outFile("atom.bin", std::ios::binary);
    if (!outFile.is_open())
    {
        std::cerr << "无法打开atom.bin" << std::endl;
        return -1;
    }
    for (int z = 0; z < zTotal; z++)
    {
        outFile.write((char*)lutStorage[z].data(), lutStorage[z].size() * sizeof(RGBAHalf));
    }
    outFile.close();
    std::cout << "3D LUT atom.bin 生成完成！" << std::endl;
    return 0;
}

int main()
{

    Write3DLUTBin();

    return 0;
}

