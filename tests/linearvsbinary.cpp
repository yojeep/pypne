#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

using namespace std;
using namespace std::chrono;

// 二分查找实现
int fsi_binary(const vector<int>& starts, int i) {
    int index = 0;
    int length = starts.size();
    while (length > 0) {
        int half = length >> 1;
        index += (starts[index + half] <= i) * (length - half);
        length = half;
    }
    return --index;
}

// 线性查找实现
int fsi_linear(const vector<int>& starts, int i) {
    int index = 0;
    while (starts[index] <= i) {
        ++index;
    }
    return --index;
}

// 生成测试数据
vector<int> generate_test_data(int n) {
    vector<int> data(n);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 100);
    
    data[0] = 0;
    for (int i = 1; i < n; ++i) {
        data[i] = data[i-1] + dis(gen);
    }
    return data;
}

// 运行性能测试
void run_benchmark(int n, int iterations = 100000) {
    vector<int> data = generate_test_data(n);
    uniform_int_distribution<> target_dis(0, data.back());
    random_device rd;
    mt19937 gen(rd());
    
    // 测试二分查找
    auto start = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        int target = target_dis(gen);
        int result = fsi_binary(data, target);
        (void)result; // 避免编译器优化掉
    }
    auto end = high_resolution_clock::now();
    auto binary_time = duration_cast<microseconds>(end - start).count();
    
    // 测试线性查找
    start = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        int target = target_dis(gen);
        int result = fsi_linear(data, target);
        (void)result; // 避免编译器优化掉
    }
    end = high_resolution_clock::now();
    auto linear_time = duration_cast<microseconds>(end - start).count();
    
    cout << "n = " << n << ":\n";
    cout << "  二分查找: " << binary_time << " μs (" << iterations << " iterations)\n";
    cout << "  线性查找: " << linear_time << " μs (" << iterations << " iterations)\n";
    cout << "  二分查找比线性查找快 " 
         << (linear_time * 100.0 / binary_time - 100) << "%\n";
    cout << "----------------------------------------\n";
}

int main() {
    // 测试不同大小的n
    vector<int> test_sizes = {5, 10, 15, 20, 25, 30, 40, 50, 75, 100};
    
    for (int n : test_sizes) {
        run_benchmark(n);
    }
    
    return 0;
}